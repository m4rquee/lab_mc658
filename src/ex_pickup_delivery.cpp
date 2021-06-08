// Project and Analysis of Algorithms
//
// Laboratorio 1
//
// Send corrections/comments to Flávio K. Miyazawa
#include "mygraphlib.h"
#include "myutils.h"
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <lemon/min_cost_arborescence.h>
#include <queue>
#include <string>

using namespace lemon;
using namespace std;

typedef vector<DNode> DNodeVector;

#define EPS 2.0

// Pickup_Delivery_Instance put all relevant information in one class.
class Pickup_Delivery_Instance {
public:
  Pickup_Delivery_Instance(Digraph &graph, DNodeStringMap &vvname,
                           DNodePosMap &posx, DNodePosMap &posy,
                           ArcValueMap &eweight, DNode &sourcenode,
                           DNode &targetnode, int &npairs, DNodeVector &pickup,
                           DNodeVector &delivery,
                           Digraph::NodeMap<DNode> &del_pickup);
  Digraph &g;
  DNodeStringMap &vname;
  DNodePosMap &px;
  DNodePosMap &py;
  ArcValueMap &weight;
  int nnodes;
  DNode &source;
  DNode &target;
  int npairs;
  DNodeVector &pickup;
  DNodeVector &delivery;
  Digraph::NodeMap<DNode> &del_pickup;
};

Pickup_Delivery_Instance::Pickup_Delivery_Instance(
    Digraph &graph, DNodeStringMap &vvname, DNodePosMap &posx,
    DNodePosMap &posy, ArcValueMap &eweight, DNode &sourcenode,
    DNode &targetnode, int &vnpairs, DNodeVector &vpickup,
    DNodeVector &vdelivery, Digraph::NodeMap<DNode> &del_pickup)
    : g(graph), vname(vvname), px(posx), py(posy), weight(eweight),
      source(sourcenode), target(targetnode), npairs(vnpairs), pickup(vpickup),
      delivery(vdelivery), del_pickup(del_pickup) {
  nnodes = countNodes(g);
}

void PrintInstanceInfo(Pickup_Delivery_Instance &P) {
  cout << endl << endl;
  cout << "Pickup Delivery Graph Informations" << endl;
  cout << "\tSource = " << P.vname[P.source] << endl;
  cout << "\tTarget = " << P.vname[P.target] << endl;
  for (int i = 0; i < P.npairs; i++) {
    cout << "\tPair pickup-->delivery: " << P.vname[P.pickup[i]] << " --> "
         << P.vname[P.delivery[i]] << endl;
  }
  cout << endl;
}

void PrintSolution(Pickup_Delivery_Instance &P, DNodeVector &Sol,
                   const string &msg) {
  // Imprime a solucao no terminal.
  cout << msg << endl << "\t";
  cout << P.vname[Sol[0]];
  for (int i = 1; i < P.nnodes; i++)
    cout << "-->" << P.vname[Sol[i]];
  cout << endl;
}

void graph_pruning(Digraph &g, const DNode &source, const DNode &target,
                   const int &npairs, const DNodeVector &pickup,
                   const DNodeVector &delivery) {
  cout << "Limpa arcos inválidos:" << endl;
  cout << "Arcos antes da limpeza: " << countArcs(g) << endl;
  // pickup[i] -> target (the last one before the target must be a delivery)
  // pickup[i] -> source (cannot go back to the source)
  for (int i = 0; i < npairs; i++)
    for (OutArcIt a(g, pickup[i]); a != INVALID; ++a)
      if (g.target(a) == target or g.target(a) == source)
        g.erase(a);
  // delivery[i] -> pickup[i] (a pickup is always before the corresponding delivery)
  // delivery[i] -> source (cannot go back to the source)
  // source -> delivery[i] (the second one after the source must be a pickup)
  for (int i = 0; i < npairs; i++) {
    for (OutArcIt a(g, delivery[i]); a != INVALID; ++a)
      if (g.target(a) == pickup[i] or g.target(a) == source)
        g.erase(a);
    for (InArcIt a(g, delivery[i]); a != INVALID; ++a)
      if (g.source(a) == source) {
        g.erase(a);
        break;
      }
  }
  // target -> v, for all nodes (target is the last node)
  for (OutArcIt a(g, target); a != INVALID; ++a)
    g.erase(a);
  // source -> target (the pickup/deliveries must be in the middle)
  for (OutArcIt a(g, source); a != INVALID; ++a)
    if (g.target(a) == target) {
      g.erase(a);
      break;
    }
  cout << "Arcos depois da limpeza: " << countArcs(g) << endl;
}

bool ReadPickupDeliveryDigraph(const string &filename, Digraph &g,
                               DNodeStringMap &vname, DNodePosMap &posx,
                               DNodePosMap &posy, ArcValueMap &weight,
                               DNode &source, DNode &target, int &npairs,
                               DNodeVector &pickup, DNodeVector &delivery,
                               Digraph::NodeMap<DNode> &del_pickup) {
  ReadDigraph(filename, g, vname, posx, posy, weight);
  int n = countNodes(g);
  DNode DN[n];
  if ((n < 4) || (n % 2)) {
    cout << "Numero de vertices " << n
         << " no grafo nao eh par ou eh menor que 4." << endl;
    return false;
  }
  npairs = (n - 2) / 2;
  pickup.resize(npairs);
  delivery.resize(npairs);
  int i = 0;
  for (DNodeIt v(g); v != INVALID; ++v) {
    DN[i] = v;
    i++;
  }

  source = DN[0];
  target = DN[1];
  for (i = 0; i < npairs; i++) {
    pickup[i] = DN[2 * i + 2];
    delivery[i] = DN[2 * i + 3];
    del_pickup[delivery[i]] = pickup[i]; // map a delivery to it's pickup
  }

  // Remove invalid arcs (7n + 2):
  graph_pruning(g, source, target, npairs, pickup, delivery);
  return true;
}

double route_cost(const Pickup_Delivery_Instance &P, const DNodeVector &Sol) {
  double cost = 0.0;
  for (int i = 1; i < P.nnodes; i++)
    for (OutArcIt a(P.g, Sol[i - 1]); a != INVALID; ++a)
      if (P.g.target(a) == Sol[i]) {
        cost += P.weight[a];
        break;
      }
  return cost;
}

// Heuristica apenas para testar a visualizacao das solucoes.
bool dummy_heuristic(Pickup_Delivery_Instance &P, int time_limit, double &LB,
                     double &UB, DNodeVector &Sol) {
  cout << "Execucao da Heuristica Boba" << endl;
  cout << "\tEsta rotina deveria respeitar o tempo de no maximo " << time_limit
       << " segundos" << endl;
  if (UB == MY_INF) { // Faz alguma coisa so' se ainda nao tem solucao
    Sol.resize(P.nnodes);
    Sol[0] = P.source; // insere o source
    for (int i = 0; i < P.npairs; i++)
      Sol[2 * i + 1] = P.pickup[i]; // insere os pickup
    for (int i = 0; i < P.npairs; i++)
      Sol[2 * i + 2] = P.delivery[i]; // insere os delivery
    Sol[2 * P.npairs + 1] = P.target; // insere o target.

    // Atualiza o UB (Upper Bound) que eh o valor da solucao
    UB = route_cost(P, Sol);
  }
  cout << endl;
  return true;
}

void _arborescence_transversal(Pickup_Delivery_Instance &P,
    MinCostArborescence<Digraph, ArcValueMap> &solver, DNodeVector &Sol,
    DNode &currNode, DCutMap &visited, map<DNode, bool> &p_visited, int pos) {
  Sol[pos] = currNode; // save the current node to the route
  if (pos == P.nnodes - 2) return; // if is in the last node before the target

  visited[currNode] = true;
  if (p_visited.count(currNode)) // mark the currNode if it's a pickup
    p_visited[currNode] = true;

  // Min arc over all neighbours:
  DNode min_next;
  double min_cost = MY_INF;
  // Min arc over all neighbours, restricted to the arborescence arcs:
  DNode min_next_arb;
  double min_cost_arb = MY_INF;

  for (OutArcIt a(P.g, currNode); a != INVALID; ++a) {
    DNode next = P.g.target(a);
    bool is_pickup = p_visited.count(next);
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
                         MinCostArborescence<Digraph, ArcValueMap> &solver) {
  map<DNode, bool> p_visited; // if each pickup has already been visited
  for (const auto &key : P.pickup) p_visited[key] = false; // init the map
  DCutMap visited(P.g, false); // if each node has already been visited

  // The target is always the last node:
  Sol[P.nnodes - 1] = P.target;
  visited[P.target] = true;

  // Transverse the graph from the source greedily guided by the arborescence:
  _arborescence_transversal(P, solver, Sol, P.source, visited, p_visited, 0);
  return route_cost(P, Sol); // Calculate the route cost
}

bool arborescence_heuristic(Pickup_Delivery_Instance &P, int time_limit, double &LB,
                            double &UB, DNodeVector &Sol) {
  // Generates the arborescence that will guide the route creation:
  MinCostArborescence<Digraph, ArcValueMap> solver(P.g, P.weight);
  solver.run(P.source); // rot the arborescence in the source
  // As a spanning digraph roted at the source this is itself a LB:
  LB = solver.arborescenceCost();

  Sol.resize(P.nnodes); // starts the solution vector
  double newUB = arborescence_transversal(P, Sol, solver);
  if (newUB < UB) { // check if this is a better solution
    UB = newUB;
    return true;
  }
  return false;
}

bool ViewPickupDeliverySolution(Pickup_Delivery_Instance &P, double &LB,
                                double &UB, DNodeVector &Sol,
                                const string &msg) {
  DigraphAttributes GA(P.g, P.vname, P.px, P.py);
  GA.SetDefaultDNodeAttrib(
      "color=LightGray style=filled width=0.2 height=0.2 fixedsize=true");
  for (ArcIt a(P.g); a != INVALID; ++a)
    GA.SetColor(a, "Invis");
  GA.SetColor(P.source, "Red"); // source and target are painted in White
  GA.SetColor(P.target, "Red");
  GA.SetShape(P.source, "star");
  GA.SetShape(P.target, "star");

  for (int i = 0; i < P.npairs; i++) { // distinguish the pickup
    GA.SetShape(P.pickup[i], "box");
  }

  if (P.npairs <= 16) { // se tiver poucos pares, dah para pintar os pares de mesma cor.
    for (int i = 0; i < P.npairs; i++) { // pinta
      GA.SetColor(P.pickup[i], ith_VisualDistinctColorName(i));
      GA.SetColor(P.delivery[i], ith_VisualDistinctColorName(i));
    }
  }
  for (int i = 1; i < P.nnodes; i++) {
    // pinta o arco Sol[i-1] -->  Sol[i]
    for (OutArcIt a(P.g, Sol[i - 1]); a != INVALID; ++a)
      if (P.g.target(a) == Sol[i]) {
        GA.SetColor(a, "Black");
        break;
      }
  }
  GA.SetLabel("Path from node " + P.vname[P.source] + " to node " +
              P.vname[P.target] + " of value " + DoubleToString(UB) +
              ". LB = " + DoubleToString(LB) + ". " + msg);
  GA.View();
  return true;
}

double *ordered_edge_prefix_sum(Pickup_Delivery_Instance &P) {
  // Init a min heap to sort the elements:
  priority_queue<double, vector<double>, greater<double>> queue;
  for (ArcIt a(P.g); a != INVALID; ++a)
    queue.push(P.weight[a]);

  // Make a prefix sum of the nnodes-1 min arcs (used as an upper bound latter):
  auto p_sum = new double[P.nnodes - 1];
  p_sum[0] = queue.top();
  queue.pop();
  for (int i = 1; i < P.nnodes - 1; i++, queue.pop())
    p_sum[i] = p_sum[i - 1] + queue.top();
  return p_sum;
}


void _exact_solution(Pickup_Delivery_Instance &P, int time_limit, double &LB,
                    double &UB, DNodeVector &Sol, double *p_sum) {
}

void exact_solution(Pickup_Delivery_Instance &P, int time_limit, double &LB,
                    double &UB, DNodeVector &Sol) {
  map<DNode, bool> p_visited; // if each pickup has already been visited
  for (const auto &key : P.pickup) p_visited[key] = false; // init the map
  DCutMap visited(P.g, false); // if each node has already been visited
  double *p_sum = ordered_edge_prefix_sum(P); // used to bound the branching
  _exact_solution(P, time_limit, LB, UB, Sol, p_sum);
}

bool Lab1(Pickup_Delivery_Instance &P, int time_limit, double &LB, double &UB,
          DNodeVector &Sol) {
  // bool improved = dummy_heuristic(P, time_limit, LB, UB, Sol);
  bool improved = arborescence_heuristic(P, time_limit, LB, UB, Sol);
  exact_solution(P, time_limit, LB, UB, Sol);
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
  DNode sourcenode, targetnode;
  Digraph::NodeMap<DNode> del_pickup(g); // map a delivery to it's pickup
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
                                 del_pickup)) {
    cout << "Erro na leitura do grafo de entrada." << endl;
    exit(EXIT_FAILURE);
  }

  Pickup_Delivery_Instance P(g, vname, px, py, weight, source, target, npairs,
                             pickup, delivery, del_pickup);
  PrintInstanceInfo(P);

  double LB = 0, UB = MY_INF; // considere MY_INF como infinito.
  DNodeVector Solucao;

  bool melhorou = Lab1(P, maxtime, LB, UB, Solucao);

  if (melhorou) {
    ViewPickupDeliverySolution(P, LB, UB, Solucao, "Solucao do Lab.");
    PrintSolution(P, Solucao, "Solucao do Lab1.");
  }
  return 0;
}
