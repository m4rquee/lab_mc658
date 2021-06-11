#ifndef LAB_MC658_PICKUP_DELIVERY_UTILS_HPP
#define LAB_MC658_PICKUP_DELIVERY_UTILS_HPP

#include "mygraphlib.h"
#include "myutils.h"
#include <chrono>
#include <queue>
#include <lemon/min_cost_arborescence.h>

#define ELAPSED ((chrono::system_clock::now() - P.start).count() / 1E9)
#define _NEW_UB_MESSAGE(SOL, MSG) {                                             \
    PrintSolution(P, SOL, MSG);                                                \
    printf("custo: %05.5f - %02.2f%% ótimo\n", UB, 100 * LB / UB);             \
  }
#define NEW_UB_MESSAGE(SOL) _NEW_UB_MESSAGE(SOL, "\nNovo UB.")

using namespace lemon;
using namespace std;

typedef chrono::time_point<chrono::system_clock> time_point;
typedef vector<DNode> DNodeVector;

// Compares two arcs based upon their weights:
struct ArcCmp : public binary_function<Arc, Arc, bool> {
  ArcValueMap &weight;

  explicit ArcCmp(ArcValueMap &weight) : weight(weight) {}

  _GLIBCXX14_CONSTEXPR
  inline bool operator()(const Arc &x, const Arc &y) const {
    return weight[x] > weight[y]; // acceding order
  }
};

typedef MinCostArborescence<Digraph, ArcValueMap> MinCostArb;
typedef priority_queue<double, vector<double>, greater<double>> min_heap;
typedef priority_queue<Arc, vector<Arc>, ArcCmp> min_arc_heap;

// Pickup_Delivery_Instance put all relevant information in one class.
class Pickup_Delivery_Instance {
public:
  Pickup_Delivery_Instance(Digraph &graph, DNodeStringMap &vvname,
                           DNodePosMap &posx, DNodePosMap &posy,
                           ArcValueMap &eweight, DNode &sourcenode,
                           DNode &targetnode, int &npairs, DNodeVector &pickup,
                           DNodeVector &delivery,
                           Digraph::NodeMap<DNode> &del_pickup,
                           DNodeBoolMap &is_pickup, int &time_limit);
  void start_counter();

  Digraph &g;
  DNodeStringMap &vname;
  DNodePosMap &px;
  DNodePosMap &py;
  ArcValueMap &weight;
  const int nnodes;
  DNode &source;
  DNode &target;
  const int npairs;
  DNodeVector &pickup;
  DNodeVector &delivery;
  Digraph::NodeMap<DNode> &del_pickup;
  map<DNode, vector<Arc>> ordered_arcs;
  DNodeBoolMap &is_pickup;
  time_point start;
  const int time_limit;
};

void PrintInstanceInfo(Pickup_Delivery_Instance &P);
void PrintSolution(Pickup_Delivery_Instance &P, DNodeVector &Sol,
                   const string &msg);

bool ReadPickupDeliveryDigraph(const string &filename, Digraph &g,
                               DNodeStringMap &vname, DNodePosMap &posx,
                               DNodePosMap &posy, ArcValueMap &weight,
                               DNode &source, DNode &target, int &npairs,
                               DNodeVector &pickup, DNodeVector &delivery,
                               Digraph::NodeMap<DNode> &del_pickup,
                               DNodeBoolMap &is_pickup);

double route_cost(const Pickup_Delivery_Instance &P, const DNodeVector &Sol);

bool ViewPickupDeliverySolution(Pickup_Delivery_Instance &P, double &LB,
                                double &UB, DNodeVector &Sol,
                                const string &msg);

#endif // LAB_MC658_PICKUP_DELIVERY_UTILS_HPP
