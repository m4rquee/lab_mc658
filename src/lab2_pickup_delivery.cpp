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
#include <string>

const long unsigned seed = 0; // seed to the random number generator
const double pe = 0.20;       // fraction of population to be the elite-set
const double pm = 0.10;       // fraction of population to be replaced by mutants
const double rhoe = 0.70;     // probability that offspring inherit an allele from elite parent
const unsigned K = 3;         // number of independent populations
const unsigned MAXT = 4;      // number of threads for parallel decoding
const unsigned X_INTVL = 100; // exchange best individuals at every 100 generations
const unsigned X_NUMBER = 2;  // exchange top 2 best
const unsigned K_MAX = 15000; // maximum value for the restart(k) strategy

inline void genArbLB(Pickup_Delivery_Instance &P, double &LB) {
  MinCostArb arb_solver(P.g, P.weight); // generates a min arborescence to derive a LB
  arb_solver.run(P.source); // root the arborescence in the source
  // As a spanning digraph rooted at the source this is itself a LB:
  LB = max(LB, arb_solver.arborescenceCost());
}

bool Lab2(Pickup_Delivery_Instance &P, double &LB, double &UB, DNodeVector &Sol) {
  bool improved = false;
  genArbLB(P, LB);

  cout << "BRKGA Variable Informations" << endl;
  const unsigned n = 2 * P.npairs; // size of chromosomes
  const unsigned p = (P.npairs > 100 ? 1 : 2) * n;   // size of population
  const unsigned k_prime = (P.npairs > 10 ? 250 : 10) * n; // candidate for k
  const unsigned k = min(k_prime, K_MAX); // restart strategy parameter
  cout << "\tsize of chromosomes       : n = " << n << endl;
  cout << "\tsize of population        : p = " << p << endl;
  cout << "\trestart strategy parameter: k = " << k << endl << endl;

  // Breaks the loop after MAX_UNCHANGED unimproved check(s):
  const unsigned MAX_UNCHANGED = P.npairs == 5 ? 1 :
                                 P.npairs == 6 ? 8 :
                                 P.npairs == 7 ? 3 :
                                 P.npairs <= 9 ? 180 : // for small instances
                                 P.npairs == 10 ? 20 : // for 10 pickups
                                 UINT_MAX;             // the rest stops by time
  const unsigned MAX_GENS =
      P.npairs >= 15 ? UINT_MAX : // on big instances limit only by time
      500000;                     // run for up to 500000 gens
  cout << "\tunchanged checks before break: MAX_UNCHANGED = " << MAX_UNCHANGED << endl;
  cout << "\tgenerations maximum          : MAX_GENS = " << MAX_GENS << endl;

  PickupDeliveryDecoder decoder(P);
  MTRand rng(seed); // initialize the random number generator
  // Initialize the BRKGA-based heuristic:
  BRKGA<PickupDeliveryDecoder, MTRand> algorithm(n, p, pe, pm, rhoe, decoder,
                                                 rng, K, MAXT);
  unsigned unchanged_checks = 0, reset_count = 0;
  DNodeVector runSol(P.nnodes); // current aux solution vector
  double best_running = MY_INF; // best fitness since last restart
  unsigned generation = 0;      // current generation
  P.start_counter(); // fixes the start time point
  do {
    algorithm.evolve(X_INTVL); // evolve the population for X_INTVL generations
    generation += X_INTVL;
    algorithm.exchangeElite(X_NUMBER); // exchange top individuals

    // Check for new UB:
    unchanged_checks++;
    reset_count++;
    double best_val_found = algorithm.getBestFitness();
    if (best_val_found < best_running) {
      reset_count = 0;
      best_running = best_val_found;
    }
    if (best_val_found < UB) {
      improved = true;
      unchanged_checks = 0;
      UB = best_val_found;
      decoder.decode(algorithm.getBestChromosome(), runSol);
      NEW_UB_MESSAGE(runSol);
      // When on large instances the local search may lead to quick improves:
      if (P.npairs >= 15)
        local_search(P, LB, UB, runSol);
    } else if (reset_count * X_INTVL >= k) { // restart(k) strategy
      cout << "\nReset depois de " << reset_count * X_INTVL << " gerações sem melhora." << endl;
      reset_count = 0;
      best_running = MY_INF;
      algorithm.reset();
    }

    int elapsed = ELAPSED;
    if (elapsed >= P.time_limit) {
      cout << "\nTempo máximo de " << P.time_limit << "s atingido." << endl;
      break;
    }
    printf("\r-> generation %d - %ds", generation, elapsed);
    cout << std::flush;
  } while (generation < MAX_GENS and unchanged_checks < MAX_UNCHANGED);

  if (unchanged_checks == MAX_UNCHANGED)
    cout << "\n\nFim após " << MAX_UNCHANGED * X_INTVL
         << " gerações sem melhora." << endl;
  else
    cout << "\n\nFim das " << MAX_GENS << " gerações." << endl;

  if (improved) Sol = runSol;
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
  DNodeVector Solucao(P.nnodes);

  bool melhorou = Lab2(P, LB, UB, Solucao);

  if (melhorou) {
    ViewPickupDeliverySolution(P, LB, UB, Solucao, "Solucao do Lab.");
    PrintSolution(P, Solucao, "\nSolucao do Lab2.");
    cout << "custo: " << UB << endl;
  }
  return 0;
}
