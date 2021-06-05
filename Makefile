
#=============================================================================
# Probably, you only need to change the attribution of the following variables:
# LEMONDIR
# GUROBI_DIR
# find the corresponding lines below. Note that for the GUROBI_DIR, there are
# two parts, depending on the operational system (Linux or OSX, see below).
# Change the paths accordingly.
#=============================================================================

#================= LEMON =====================================================
# if you have the lemon package installed in a specific folder, change the
# following line:
LEMONDIR = ./lemon/lemon-1.3.1

LEMONINC  = -I$(LEMONDIR)/include
LEMONLIB  = -L$(LEMONDIR)/lib   -lemon 
#================= GUROBI =====================================================
# If Gurobi is installed, there is a program called gurobi_cl that can be used to
# detect the version
VERSION := $(shell gurobi_cl --version | cut -c 26,28,30 | head -n 1)
FLAGVERSION := $(shell gurobi_cl --version | cut -c 26,28 | head -n 1)

HOMEDIR = .

# if your operational system is OSX (e.g., for mac computers)
ifeq ($(shell uname), Darwin)
        $(info )
        $(info *)
        $(info *        Makefile for MAC OS environment)
        $(info *)
	PLATFORM = mac64
	ifneq ($(RELEASE), 11)
                # The next variable must be used to compile *and* link the obj. codes
		CPPSTDLIB = -stdlib=libc++
	else
                $(info )
                $(info *    Gurobi library is not compatible with code)
                $(info *    generated by clang c++11 in MAC OS.)
                $(info *    Please, verify if it is compatible now.)
                $(info *)
                $(info *    >>>>> Aborted <<<<<)
                $(info *)
                $(error *)
		CPPSTDLIB = -stdlib=libc++ -std=c++11
	endif
	CC      = g++
	#CC_ARGS    = -Wall -m64 -O3 -Wall $(CPPSTDLIB)  -Wno-c++11-extensions -Wunused-local-typedefs
	#CC_ARGS    = -m64 -g -ferror-limit=2 -Wall -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0  -Wunused-local-typedefs
	CC_ARGS    = -m64 -O2 -ferror-limit=2 -Wall -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 -Wunused-local-typedefs
	RELEASE := $(shell uname -r | cut -f 1 -d .)
	CC_LIB   = -lm -lpthread $(CPPSTDLIB)
	GUROBI_DIR = /Library/gurobi$(VERSION)/$(PLATFORM)

	# if your operational system is Linux
else
        $(info )
        $(info *)
        $(info *        Makefile for LINUX environment)
        $(info *)
	PLATFORM = linux64
	CC      = g++
	#CC_ARGS    = -m64 -O2 -Wall -std=c++11 -Wunused-local-typedefs
	CC_ARGS    = -m64 -O2 -Wall -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0  -Wunused-local-typedefs
	RELEASE := $(shell uname -r | cut -f 1 -d .)
	CC_LIB   = -lm -lpthread
	GUROBI_DIR = /opt/gurobi$(VERSION)/$(PLATFORM)
endif
GUROBI_INC = -I$(GUROBI_DIR)/include
GUROBI_LIB = -L$(GUROBI_DIR)/lib  -lgurobi_c++ -lgurobi$(FLAGVERSION)  $(CPPSTDLIB)
#===============================================================================


HOMEDIR_INC = $(HOMEDIR)/include
HOMEDIR_LIB = $(HOMEDIR)/lib
HOMEDIR_SRC = $(HOMEDIR)/src
HOMEDIR_OBJ = $(HOMEDIR)/obj
HOMEDIR_BIN = $(HOMEDIR)/bin

#---------------------------------------------
# define includes and libraries
#INC = $(GUROBI_INC)  $(LEMONINC) -I$(HOMEDIR_INC) 
#LIB = $(CC_LIB) $(GUROBI_LIB) $(LEMONLIB) -L$(HOMEDIR_LIB) 
INC = $(LEMONINC) -I$(HOMEDIR_INC) 
LIB = $(CC_LIB) $(LEMONLIB) -L$(HOMEDIR_LIB) 

_MLS = mygraphlib.cpp geompack.cpp myutils.cpp mycolor.cpp deprecated.cpp
_MLO = $(_MLS:.cpp=.o)
MYLIB_SRC = $(patsubst %,$(HOMEDIR_SRC)/%,$(_MLS))
MYLIB_OBJ = $(patsubst %,$(HOMEDIR_OBJ)/%,$(_MLO))

#_EX = ex_shortestpath.cpp ex_steiner_undirected.cpp ex_greedycoloring.cpp ex_graphtable.cpp ex_mycolor.cpp ex_coloring.cpp ex_visualcolors.cpp ex_cvrp.cpp ex_tsp_mtz.cpp test_tsp.cpp ex_basics_graph.cpp ex_fractional_packing.cpp ex_knapsack.cpp ex_tsp.cpp ex_fractional_tsp.cpp  generate_cflp.cpp generate_random_euclidean_graph.cpp generate_random_graph.cpp generate_bipartite_vertex_cover.cpp generate_triangulated_digraph.cpp generate_triangulated_graph.cpp ex_steiner_directed.cpp generate_steiner_file.cpp ex_kpaths.cpp ex_cflp.cpp ex_bipartite_matching.cpp ex_bipartite_matching2.cpp ex_perfect_matching_general_graphs.cpp ex_viewgraph.cpp ex_viewdigraph.cpp ex_vertex_cover_in_bipartite_graph.cpp 
_EX = ex_shortestpath.cpp ex_greedycoloring.cpp ex_graphtable.cpp ex_mycolor.cpp ex_visualcolors.cpp ex_basics_graph.cpp generate_cflp.cpp generate_random_euclidean_graph.cpp generate_random_graph.cpp generate_bipartite_vertex_cover.cpp generate_triangulated_digraph.cpp generate_triangulated_graph.cpp ex_viewgraph.cpp ex_viewdigraph.cpp  generate_pickup_delivery_digraph.cpp ex_pickup_delivery.cpp
_OB = $(_EX:.cpp=.o)
_BN = $(_EX:.cpp=.e)

# complete the EX's names with the path
EXAMPLES_SRC = $(patsubst %,$(HOMEDIR_SRC)/%,$(_EX))
EXAMPLES_OBJ = $(patsubst %,$(HOMEDIR_OBJ)/%,$(_OB))
EXAMPLES_BIN = $(patsubst %,$(HOMEDIR_BIN)/%,$(_BN))

all: $(EXAMPLES_OBJ) $(EXAMPLES_BIN) $(EXAMPLES_SRC) $(MYLIB_SRC) $(MYLIB_OBJ) $(HOMEDIR_LIB)/mylib.a


$(HOMEDIR_LIB)/mylib.a: $(MYLIB_SRC) $(MYLIB_OBJ)
	#libtool -o $@ $(MYLIB_OBJ)
	ar cru $@ $(MYLIB_OBJ)

$(HOMEDIR_BIN)/%.e: $(HOMEDIR_OBJ)/%.o $(HOMEDIR_LIB)/mylib.a
	$(CC) $(CC_ARGS) $^ -o $@ $(LIB) $(INC) 

$(HOMEDIR_OBJ)/%.o: $(HOMEDIR_SRC)/%.cpp
	$(CC) $(CC_ARGS) -c $^ -o $@ $(INC) 



clean:
	rm -f $(HOMEDIR)/*~ $(HOMEDIR_BIN)/*.e $(HOMEDIR_OBJ)/*.o $(HOMEDIR_LIB)/*.a $(HOMEDIR_SRC)/*~ $(HOMEDIR_INC)/*~ 
