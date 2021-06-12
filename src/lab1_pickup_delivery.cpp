// Project and Analysis of Algorithms
//
// Laboratorio 1
//
// Send corrections/comments to Flávio K. Miyazawa
#include "pickup_delivery_utils.hpp"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <string>

#define EPS 2.0

void _arborescence_transversal(Pickup_Delivery_Instance &P, MinCostArb &solver,
                               DNodeVector &Sol, DNode &currNode,
                               DNodeBoolMap &visited,
                               map<DNode, bool> &p_visited, int pos) {
  Sol[pos] = currNode;             // save the current node to the route
  if (pos == P.nnodes - 2) return; // if is in the last node before the target

  // Mark as visited:
  visited[currNode] = true;
  if (P.is_pickup[currNode]) // mark the currNode if it's a pickup
    p_visited[currNode] = true;

  // Min arc over all neighbours:
  DNode min_next;
  double min_cost = MY_INF;
  // Min arc over all neighbours, restricted to the arborescence arcs:
  DNode min_next_arb;
  double min_cost_arb = MY_INF;

  for (const Arc &a : P.ordered_arcs[currNode]) {
    DNode next = P.g.target(a);
    bool is_pickup = P.is_pickup[next];
    // Can only go to next if not visited, and it's a pickup or it's
    // corresponding pickup has already been visited:
    if (!visited[next] && (is_pickup || p_visited[P.del_pickup[next]])) {
      if (P.weight[a] < min_cost) { // min arc over all neighbours
        min_next = next;
        min_cost = P.weight[a];
      }
      if (P.weight[a] < min_cost_arb && solver.arborescence(a)) { // restricted min
        min_next_arb = next;
        min_cost_arb = P.weight[a];
        break; // the arcs are sorted, so the first one is the minimum
      }
    }
  }

  // Choose the min arc giving preference to the restricted one:
  DNode &nextNode = min_next;
  if (min_cost_arb != MY_INF && min_cost_arb <= min_cost * EPS)
    nextNode = min_next_arb;
  _arborescence_transversal(P, solver, Sol, nextNode, visited, p_visited, ++pos);
}

double arborescence_transversal(Pickup_Delivery_Instance &P, DNodeVector &Sol,
                                MinCostArb &solver) {
  map<DNode, bool> p_visited; // if each pickup has already been visited
  for (const auto &key : P.pickup) p_visited[key] = false; // init the map
  DNodeBoolMap visited(P.g, false); // if each node has already been visited

  // The target is always the last node:
  Sol[P.nnodes - 1] = P.target;
  visited[P.target] = true;

  // Transverse the graph from the source greedily guided by the arborescence:
  _arborescence_transversal(P, solver, Sol, P.source, visited, p_visited, 0);
  return route_cost(P, Sol); // Calculate the route cost
}

bool arborescence_heuristic(Pickup_Delivery_Instance &P, double &LB, double &UB,
                            DNodeVector &Sol, MinCostArb &solver) {
  DNodeVector arbSol(P.nnodes); // starts the solution vector
  double newUB = arborescence_transversal(P, arbSol, solver);
  if (newUB < UB) { // check if this is a better solution
    UB = newUB;
    Sol = arbSol;
    _NEW_UB_MESSAGE(Sol, "Novo UB pela heuristica da arborescência.");
    cout << "LB: " << LB << endl << endl;
    return true;
  }
  return false;
}

double *ordered_edge_prefix_sum(Pickup_Delivery_Instance &P, MinCostArb &solver) {
  // Init a min heap to sort the elements:
  min_heap queue;
  for (ArcIt a(P.g); a != INVALID; ++a)
    if (solver.arborescence(a)) // restricted arcs to increase the bounds
      queue.push(P.weight[a]);

  // Make a prefix sum of the nnodes-1 min arcs (used as an upper bound latter):
  auto p_sum = new double[P.nnodes - 1];
  p_sum[0] = queue.top();
  queue.pop();
  for (int i = 1; i < P.nnodes - 1; i++, queue.pop())
    p_sum[i] = p_sum[i - 1] + queue.top();
  return p_sum;
}

bool _exact_solution(Pickup_Delivery_Instance &P, double &LB, double &UB,
                     DNodeVector &currSol, DNodeVector &bestSol, double *p_sum,
                     DNode &currNode, DNodeBoolMap &visited,
                     map<DNode, bool> &p_visited, int pos, double curr_weight,
                     double &bound, bool &broke) {
  bool improved = false;
  int elapsed = ELAPSED;
  // Debug message:
  if (P.nnodes - pos >= 20) { // cut when the frequency speeds up
    printf("\r-> nó %4s - nível %d - %ds", P.vname[currNode].data(), pos, elapsed);
    cout << std::flush;
  }
  if (elapsed >= P.time_limit) {
    broke = true;
    return improved;
  }

  currSol[pos] = currNode; // save the current node to the route
  // Evaluate the newly created route when it's done:
  if (pos == P.nnodes - 2) {
    // Add the last arc weight:
    for (OutArcIt a(P.g, currNode); a != INVALID; ++a)
      if (P.g.target(a) == P.target) {
        curr_weight += P.weight[a];
        break;
      }
    // Test if it's the new best solution:
    if (curr_weight < UB) {
      UB = curr_weight;
      bestSol = currSol;
      NEW_UB_MESSAGE(bestSol);
      improved = true;
      // When on large instances the local search may lead to quick improves:
      if (P.npairs >= 15)
        local_search(P, LB, UB, bestSol);
    }
    if (UB < bound) bound = UB; // can lower the bound even more
    return improved;
  }

  for (const Arc &a : P.ordered_arcs[currNode]) {
    DNode next = P.g.target(a);
    bool is_pickup = P.is_pickup[next];
    // Can only go to next if not visited, and it's a pickup or it's
    // corresponding pickup has already been visited:
    if (!visited[next] && (is_pickup || p_visited[P.del_pickup[next]])) {
      visited[next] = true;
      if (is_pickup) p_visited[next] = true;

      int remaining_arcs = P.nnodes - pos - 2; // used pos + 1 arcs already
      double min_weight_for_remaining_arcs = p_sum[remaining_arcs - 1];
      double new_weight = curr_weight + P.weight[a];
      // Bound:
      if (new_weight < bound - min_weight_for_remaining_arcs) // branch only if can beat the bound
        // Branch:
        improved |= _exact_solution(P, LB, UB, currSol, bestSol, p_sum, next, visited,
                            p_visited, pos + 1, new_weight, bound, broke);

      if (is_pickup) p_visited[next] = false;
      visited[next] = false;
      if (LB == UB) return true; // found an optimal solution
    }
  }
  return improved;
}

bool exact_solution(Pickup_Delivery_Instance &P, double &LB, double &UB,
                    DNodeVector &Sol, MinCostArb &solver) {
  map<DNode, bool> p_visited; // if each pickup has already been visited
  for (const auto &key : P.pickup) p_visited[key] = false; // init the map
  DNodeBoolMap visited(P.g, false); // if each node has already been visited
  double *p_sum = ordered_edge_prefix_sum(P, solver); // used to bound the branching

  // The target is always the last node:
  DNodeVector currSol(P.nnodes);
  currSol[P.nnodes - 1] = P.target;

  // Mark the source/target as visited to block visiting them in the search:
  visited[P.source] = true;
  visited[P.target] = true;

  bool improved = false, broke = false;
  double bound = 1.2 * LB; // used as a fake UB in order to possibly raise the LB
  // If the lower bound is too low tries to raise it (only if is fast):
  if (bound < UB and P.npairs <= 15) {
    cout << "-----> Tenta melhorar o LB:" << endl;
    improved |= _exact_solution(P, LB, UB, currSol, Sol, p_sum, P.source,
                                visited, p_visited, 0, 0, bound, broke);
    cout << "\r";
    if (UB > bound && !broke) { // no solution could beat this bound
      LB = bound; // if no solution was found then all valid solutions have cost >= bound
      cout << "Novo LB - " << LB << endl;
    }
  }
  if (LB != UB) { // actually solves optimally
    cout << endl << "-----> Resolve otimamente:" << endl;
    improved |= _exact_solution(P, LB, UB, currSol, Sol, p_sum, P.source,
                                visited, p_visited, 0, 0, UB, broke);
    cout << endl << endl;
  }
  if (!broke)
    LB = UB; // After testing all possibilities the Sol must be optimal
  else
    cout << endl << "Tempo máximo de " << P.time_limit << "s atingido." << endl;
  return improved;
}

bool Lab1(Pickup_Delivery_Instance &P, double &LB, double &UB, DNodeVector &Sol) {
  P.start_counter(); // fixes the start time point

  // Generates the arborescence that will guide the route creation:
  MinCostArb arb_solver(P.g, P.weight);
  arb_solver.run(P.source); // root the arborescence in the source
  // As a spanning digraph rooted at the source this is itself a LB:
  LB = max(LB, arb_solver.arborescenceCost());

  bool improved = arborescence_heuristic(P, LB, UB, Sol, arb_solver);

  if (LB != UB) // if can improve
    improved |= local_search(P, LB, UB, Sol);

  if (LB != UB) // if can improve
    improved |= exact_solution(P, LB, UB, Sol, arb_solver);
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
  int seed = 0;
  srand48(seed);

  // uncomment one of these lines to change default pdf reader, or insert new
  // one set_pdfreader("open");    // pdf reader for Mac OS X
  // set_pdfreader("xpdf");    // pdf reader for Linux
  // set_pdfreader("evince");  // pdf reader for Linux
  // set_pdfreader("open -a Skim.app");
  set_pdfreader("xdg-open"); // the Linux will choose the default one
  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc != 3) {
    cout << endl
         << "Laboratorio de MC658: Rota com coleta e entrega de peso minimmo,"
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

  double LB = 0, UB = MY_INF; // considere MY_INF como infinito.
  DNodeVector Solucao;

  bool melhorou = Lab1(P, LB, UB, Solucao);

  if (melhorou) {
    ViewPickupDeliverySolution(P, LB, UB, Solucao, "Solucao do Lab.");
    PrintSolution(P, Solucao, "Solucao do Lab1.");
    cout << "custo: " << UB << endl;
  }
  return 0;
}
