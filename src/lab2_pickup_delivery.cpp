// Project and Analysis of Algorithms
//
// Laboratorio 2
//
// Send corrections/comments to Flávio K. Miyazawa
#include "PickupDeliveryDecoder.hpp"
#include "pickup_delivery_utils.hpp"
#include <BRKGA.h>
#include <MTRand.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/preflow.h>
#include <string>

const long unsigned seed = 0; // seed to the random number generator
const unsigned p = 1000;      // size of population
const double pe = 0.20;       // fraction of population to be the elite-set
const double pm = 0.10;       // fraction of population to be replaced by mutants
const double rhoe = 0.70;     // probability that offspring inherit an allele from elite parent
const unsigned K = 3;         // number of independent populations
const unsigned MAXT = 4;      // number of threads for parallel decoding

const unsigned MAX_GENS = 1000;	// run for 1000 gens
const unsigned X_INTVL = 100;	// exchange best individuals at every 100 generations
const unsigned X_NUMBER = 2;	// exchange top 2 best

bool Lab2(Pickup_Delivery_Instance &P, double &LB, double &UB, DNodeVector &Sol) {
  bool improved = false;
  P.start_counter(); // fixes the start time point

  const unsigned n = 2 * P.npairs; // size of chromosomes
  PickupDeliveryDecoder decoder;
  MTRand rng(seed); // initialize the random number generator
  // Initialize the BRKGA-based heuristic:
  BRKGA<PickupDeliveryDecoder, MTRand> algorithm(n, p, pe, pm, rhoe, decoder,
                                                 rng, K, MAXT);
  unsigned generation = 0; // current generation
  do {
    algorithm.evolve(); // evolve the population for one generation

    if (++generation % X_INTVL == 0)
      algorithm.exchangeElite(X_NUMBER); // exchange top individuals
  } while (generation < MAX_GENS);

  double best_val_found = algorithm.getBestFitness();
  if (best_val_found < UB) {
    UB = best_val_found;
    decoder.decode(algorithm.getBestChromosome(), Sol);
    improved = true;
  }
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

  bool melhorou = Lab2(P, LB, UB, Solucao);

  if (melhorou) {
    ViewPickupDeliverySolution(P, LB, UB, Solucao, "Solucao do Lab.");
    PrintSolution(P, Solucao, "Solucao do Lab2.");
    cout << "custo: " << UB << endl;
  }
  return 0;
}
