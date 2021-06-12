#include "pickup_delivery_utils.hpp"

Pickup_Delivery_Instance::Pickup_Delivery_Instance(
    Digraph &graph, DNodeStringMap &vvname, DNodePosMap &posx,
    DNodePosMap &posy, ArcValueMap &eweight, DNode &sourcenode,
    DNode &targetnode, int &vnpairs, DNodeVector &vpickup,
    DNodeVector &vdelivery, Digraph::NodeMap<DNode> &del_pickup,
    DNodeBoolMap &is_pickup, int &time_limit)
    : g(graph), vname(vvname), px(posx), py(posy), weight(eweight),
      nnodes(2 * vnpairs + 2), source(sourcenode), target(targetnode),
      npairs(vnpairs), pickup(vpickup), delivery(vdelivery),
      del_pickup(del_pickup), is_pickup(is_pickup), time_limit(time_limit) {
  // Store the out arcs of each node sorted by weight:
  ArcCmp arcCmp(weight);           // arc comparator based on this weight map
  min_arc_heap sorting_heap(arcCmp); // aux heap used for sorting
  // Sort the nodes out arcs for each node:
  for (DNodeIt n(g); n != INVALID; ++n) {
    ordered_arcs[n].reserve(npairs);
    for (OutArcIt a(g, n); a != INVALID; ++a) // add all out arcs to the heap
      sorting_heap.push(a);
    while (!sorting_heap.empty()) { // sort the arcs by popping the heap
      ordered_arcs[n].push_back(sorting_heap.top());
      sorting_heap.pop();
    }
  }
}

void Pickup_Delivery_Instance::start_counter() {
  start = chrono::system_clock::now();
}

void PrintInstanceInfo(Pickup_Delivery_Instance &P) {
  cout << endl << endl;
  cout << "Pickup Delivery Graph Informations" << endl;
  cout << "\tTime limit = " << P.time_limit << "s" << endl;
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
                               Digraph::NodeMap<DNode> &del_pickup,
                               DNodeBoolMap &is_pickup) {
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
    is_pickup[pickup[i]] = true;
    del_pickup[delivery[i]] = pickup[i]; // map a delivery to it's pickup
  }

  // Remove invalid arcs (7n + 2):
  graph_pruning(g, source, target, npairs, pickup, delivery);
  return true;
}

double route_cost(const Pickup_Delivery_Instance &P, const DNodeVector &Sol) {
  double cost = 0.0;
  bool valid = false;
  for (int i = 1; i < P.nnodes; i++, valid = false) {
    for (OutArcIt a(P.g, Sol[i - 1]); a != INVALID; ++a)
      if (P.g.target(a) == Sol[i]) {
        cost += P.weight[a];
        valid = true;
        break;
      }
    if (!valid) return MY_INF; // the route is invalid
  }
  return cost;
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

bool can_swap(Pickup_Delivery_Instance &P, DNodeVector &Sol, int i, int j) {
  bool i_is_pickup = P.is_pickup[Sol[i]];
  bool j_is_pickup = P.is_pickup[Sol[j]];
  if (i_is_pickup and j_is_pickup) {
    for (int k = i + 1; k < j; k++)
      if (!P.is_pickup[Sol[k]] and Sol[i] == P.del_pickup[Sol[k]])
        return false; // i is delivered at k
  } else if (i_is_pickup and !j_is_pickup) {
    if (Sol[i] == P.del_pickup[Sol[j]]) return false; // i is delivered at j
    for (int k = i + 1; k < j; k++)
      if ((!P.is_pickup[Sol[k]] and Sol[i] == P.del_pickup[Sol[k]]) or
          Sol[k] == P.del_pickup[Sol[j]])
        return false; // i is delivered at k or k is delivered at j
  } else if (!i_is_pickup and j_is_pickup)
    return true; // no problem swapping
  else // (!i_is_pickup and !j_is_pickup):
    for (int k = i + 1; k < j; k++)
      if (Sol[k] == P.del_pickup[Sol[j]])
        return false; // k is delivered at j
  return true; // can swap the nodes
}

bool _local_search(Pickup_Delivery_Instance &P, double &LB, double &UB,
                   DNodeVector &Sol) {
  bool improved = false;
  double new_cost;
  int n = P.nnodes - 2;
  // Search from end to begin because the heavier arcs are there:
  for (int i = n; i >= 1; i--)
    for (int j = i - 1; j >= 1; j--)
      if (can_swap(P, Sol, j, i)) { // the first index must be the smaller one
        swap(Sol[i], Sol[j]);
        new_cost = route_cost(P, Sol);
        if (new_cost < UB) { // found a better solution
          improved = true;
          UB = new_cost;
          NEW_UB_MESSAGE(Sol);
        } else
          swap(Sol[i], Sol[j]); // reset
      }
  return improved;
}

bool local_search(Pickup_Delivery_Instance &P, double &LB, double &UB,
                  DNodeVector &Sol) {
  bool improved = false, aux;
  cout << "-----> Fazendo uma busca local." << endl;
  while ((aux = _local_search(P, LB, UB, Sol))) {
    improved |= aux;
    int elapsed = ELAPSED;
    if (elapsed >= P.time_limit) {
      cout << endl << "Tempo máximo de " << P.time_limit << "s atingido." << endl;
      break;
    }
  }
  cout << endl;
  return improved;
}
