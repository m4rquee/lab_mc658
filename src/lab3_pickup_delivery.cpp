// Project and Analysis of Algorithms
//
// Laboratorio 3
//
// Send corrections/comments to Flávio K. Miyazawa
#include "pickup_delivery_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <gurobi_c++.h>
#include <iostream>
#include <lemon/list_graph.h>
#include <solver.h>
#include <string>

int LAZY_ADD = INT32_MAX; // saves memory by limiting the calls to addLazy
const long unsigned seed = 42; // seed to the random number generator

class SubCycleElim : public GRBCallback {
  const unsigned &M;
  Pickup_Delivery_Instance &P;
  Digraph::ArcMap<GRBVar> &x_e;
  Digraph::NodeMap<GRBVar> &p_v;
  double (GRBCallback::*solution_value)(GRBVar) = nullptr;

public:
  SubCycleElim(Pickup_Delivery_Instance &p, Digraph::ArcMap<GRBVar> &x_e,
               Digraph::NodeMap<GRBVar> &p_v, const unsigned int &m)
      : M(m), P(p), x_e(x_e), p_v(p_v) {}

protected:
  void callback() override {
    // -------------------------------------------------------------------------
    // Get the correct function to obtain the values of the lp variables:
    if (where == GRB_CB_MIPSOL) // if this condition is true, all variables are integer
      solution_value = &SubCycleElim::getSolution;
    else if (where == GRB_CB_MIPNODE && // node with optimal fractional solution
             getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)
      solution_value = &SubCycleElim::getNodeRel;
    else
      return; // return, as this code do not take advantage of the other options

    try {
      int constrCount = 0;
      for (ArcIt e(P.g); e != INVALID; ++e) {
        const DNode &u = P.g.source(e), &v = P.g.target(e);
        const int u_pos = (int)std::round((this->*solution_value)(p_v[u])),
                  v_pos = (int)std::round((this->*solution_value)(p_v[v]));
        if ((this->*solution_value)(x_e[e]) > 1 - MY_EPS) // this edge is used
          if (abs(v_pos - u_pos) > 1 - MY_EPS) { // arc between nonadjacent nodes
            // Arcs only between adjacent nodes:
            addLazy(p_v[v] - p_v[u] + M * (1 - x_e[e]) >= x_e[e]);
            addLazy(p_v[v] - p_v[u] - M * (1 - x_e[e]) <= x_e[e]);
            constrCount += 2;
            if (constrCount >= LAZY_ADD) break; // adds up to LAZY_ADD restrictions
          }
      }
    } catch (std::exception &e) {
      cout << "Error during callback: " << e.what() << endl;
    }
  }
};

inline bool arb_heuristic(Pickup_Delivery_Instance &P, double &LB, double &UB, DNodeVector &Sol) {
  // Generates the arborescence that will guide the route creation:
  MinCostArb arb_solver(P.g, P.weight);
  arb_solver.run(P.source); // root the arborescence in the source
  // As a spanning digraph rooted at the source this is itself a LB:
  LB = max(LB, arb_solver.arborescenceCost());
  bool improved = arborescence_heuristic(P, LB, UB, Sol, arb_solver);
  return local_search(P, LB, UB, Sol) or improved;
}

bool Lab3(Pickup_Delivery_Instance &P, double &LB, double &UB, DNodeVector &Sol) {
  P.start_counter();
  bool improved = arb_heuristic(P, LB, UB, Sol);

  // Gurobi ILP problem setup:
  GRBEnv env = GRBEnv();
  env.set(GRB_IntParam_Seed, seed);
  env.set(GRB_DoubleParam_TimeLimit, P.time_limit);
  env.set(GRB_DoubleParam_Cutoff, UB); // set the know UB
  env.set(GRB_IntParam_LazyConstraints, 1); // enable lazy constraints
  GRBModel model = GRBModel(env);
  model.set(GRB_StringAttr_ModelName, "Pickup Delivery Route");
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

  // ILP solver parameters: ----------------------------------------------------
  if (P.npairs >= 100) { // reduce memory usage
    LAZY_ADD = 100;
    model.set(GRB_IntParam_Threads, 1);
  }

  if (P.npairs >= 20) { // focus only on new UBs
    model.set(GRB_IntParam_MIPFocus, GRB_MIPFOCUS_FEASIBILITY);
    model.set(GRB_IntParam_Cuts, GRB_CUTS_AGGRESSIVE);
    model.set(GRB_IntParam_Presolve, GRB_PRESOLVE_AGGRESSIVE);
    model.set(GRB_IntParam_MinRelNodes, 500);
    model.set(GRB_IntParam_PumpPasses, 5000);
    model.set(GRB_IntParam_ZeroObjNodes, 500);
    model.set(GRB_DoubleParam_Heuristics, 0.3);
  }

  // ILP problem variables: ----------------------------------------------------
  Digraph::ArcMap<GRBVar> x_e(P.g); // binary variables for each arc
  // Binary variables that indicates if a node is in a position:
  Digraph::NodeMap<GRBVar> p_v(P.g); // position of each node in the solution
  for (ArcIt e(P.g); e != INVALID; ++e) {
    char name[100];
    sprintf(name, "x_(%s,%s)", P.vname[P.g.source(e)].c_str(),
            P.vname[P.g.target(e)].c_str());
    x_e[e] = model.addVar(0.0, 1.0, P.weight[e], GRB_BINARY, name);
  }
  for (DNodeIt v(P.g); v != INVALID; ++v) {
    char name[100];
    sprintf(name, "p_%s", P.vname[v].c_str());
    p_v[v] = model.addVar(0.0, P.nnodes - 1, 0, GRB_INTEGER, name);
  }
  model.update(); // run update to use model inserted variables

  // ILP problem restrictions: -------------------------------------------------
  cout << "Adding the model restrictions:" << endl;
  // The source and the target are fixed:
  model.addConstr(p_v[P.source] == 0);
  model.addConstr(p_v[P.target] == P.nnodes - 1);

  int constrCount = 0;
  for (const auto &delivery : P.delivery) {
    model.addConstr(p_v[P.del_pickup[delivery]] <= p_v[delivery] - 1); // the ith pickup shows up before the ith delivery
    constrCount++;
  }
  cout << "-> the ith pickup shows up before the ith delivery - " << constrCount << " constrs" << endl;

  GRBLinExpr s_out_degree_expr;
  for (OutArcIt e(P.g, P.source); e != INVALID; ++e) s_out_degree_expr += x_e[e];
  model.addConstr(s_out_degree_expr == 1); // the source is the first node
  GRBLinExpr t_in_degree_expr;
  for (InArcIt e(P.g, P.target); e != INVALID; ++e) t_in_degree_expr += x_e[e];
  model.addConstr(t_in_degree_expr == 1); // the target is the last node

  constrCount = 0;
  for (DNodeIt v(P.g); v != INVALID; ++v) {
    if (v == P.source or v == P.target) continue;
    GRBLinExpr out_degree_expr, in_degree_expr;
    for (OutArcIt e(P.g, v); e != INVALID; ++e) out_degree_expr += x_e[e];
    for (InArcIt e(P.g, v); e != INVALID; ++e) in_degree_expr += x_e[e];
    // The in/out-degree of each internal node is one:
    model.addConstr(out_degree_expr == 1);
    model.addConstr(in_degree_expr == 1);
    constrCount += 2;
  }
  cout << "-> the in/out-degree of each internal node is one - " << constrCount << " constrs" << endl;

  unsigned M = P.nnodes;
  SubCycleElim cb(P, x_e, p_v, M);
  model.setCallback(&cb);
  cout << "-> arcs only between adjacent nodes - adding a callback" << endl;
  cout << "-> other - " << 4 << " constrs" << endl;

  // ILP solving: --------------------------------------------------------------
  model.optimize(); // trys to solve optimally within the time limit
  LB = max(LB, model.get(GRB_DoubleAttr_ObjBound)); // updates the LB
  if (model.get(GRB_IntAttr_SolCount) > 0) {  // a better solution was found
    improved = true;
    UB = GetModelValue(model);
    for (DNodeIt v(P.g); v != INVALID; ++v) // saves this better solution
      Sol[(unsigned)std::round(p_v[v].get(GRB_DoubleAttr_X))] = v;
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) // solved optimally
      LB = UB;
    else { // can still improve
      NEW_UB_MESSAGE(Sol); // print the ILP solution
      cout << endl;
      local_search(P, LB, UB, Sol);
    }
    NEW_UB_MESSAGE(Sol);
  }
  cout << "Novo LB - " << LB << endl;

  return improved;
}

int main(int argc, char *argv[]) {
  int maxtime;
  Digraph g; // graph declaration
  string digraph_filename, source_node_name, target_node_name;
  DNodeStringMap vname(g);  // name of graph nodes
  DNodePosMap px(g), py(g); // xy-coodinates for each node
  DNodeColorMap vcolor(g);  // color of nodes
  ArcStringMap aname(g);    // name for graph arcs
  ArcColorMap ecolor(g);    // color of edges
  ArcValueMap lpvar(g);     // used to obtain the contents of the LP variables
  ArcValueMap weight(g);    // edge weights
  vector<DNode> V;
  Digraph::NodeMap<DNode> del_pickup(g); // map a delivery to it's pickup
  DNodeBoolMap is_pickup(g, false); // used to quickly check if a node is a pickup
  srand48(seed);

  // uncomment one of these lines to change default pdf reader, or insert new
  // one set_pdfreader("open");    // pdf reader for Mac OS X
  // set_pdfreader("xpdf");    // pdf reader for Linux
  // set_pdfreader("evince");  // pdf reader for Linux
  // set_pdfreader("open -a Skim.app");
  set_pdfreader("xdg-open"); // the Linux will choose the default one
  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc < 3) {
    cout << endl
         << "Laboratorio de MC658: Rota com coleta e entrega de peso minimo,"
         << endl
         << "the st-shortest path problem." << endl
         << endl
         << "Usage: " << argv[0]
         << "  <pickup_delivery_digraph_filename> <maximum_time_sec>" << endl
         << endl;
    cout << "Example:" << endl
         << "\t" << argv[0] << " "
         << getpath(argv[0]) + "../instances/pickup_delivery_5.dig 10" << endl
         << endl
         << "\t" << argv[0] << " "
         << getpath(argv[0]) + "../instances/pickup_delivery_10.dig 100" << endl
         << endl;
    exit(0);
  }

  digraph_filename = argv[1];
  maxtime = atoi(argv[2]);
  double LB = 0, UB = MY_INF; // considere MY_INF como infinito.
  if (argc >= 4)
    LB = atof(argv[3]);
  if (argc >= 5)
    UB = atof(argv[4]);
  DNodeVector pickup, delivery;
  DNode source, target;
  int npairs;

  if (!ReadPickupDeliveryDigraph(digraph_filename, g, vname, px, py, weight,
                                 source, target, npairs, pickup, delivery,
                                 del_pickup, is_pickup)) {
    cout << "Erro na leitura do grafo de entrada." << endl;
    exit(EXIT_FAILURE);
  }

  Pickup_Delivery_Instance P(g, vname, px, py, weight, source, target, npairs,
                             pickup, delivery, del_pickup, is_pickup, maxtime);
  PrintInstanceInfo(P);

  DNodeVector Solucao(P.nnodes);

  try {
    bool melhorou = Lab3(P, LB, UB, Solucao);

    if (melhorou) {
      ViewPickupDeliverySolution(P, LB, UB, Solucao, "Solucao do Lab.");
      PrintSolution(P, Solucao, "\nSolucao do Lab3.");
      cout << "custo: " << UB << endl;
    }
  } catch (std::exception &e) {
    cerr << "\nException: " << e.what() << endl;
    return 1;
  }
  return 0;
}
