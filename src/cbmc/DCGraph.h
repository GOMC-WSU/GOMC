/*******************************************************************************
GPU OPTIMIZED MONTE CARLO (GOMC) 1.0 (Serial version)
Copyright (C) 2015  GOMC Group
A copy of the GNU General Public License can be found in the COPYRIGHT.txt
along with this program, also can be found at <http://www.gnu.org/licenses/>.
********************************************************************************/
#ifndef DCGRAPH_H
#define DCGRAPH_H
#include "../CBMC.h"
#include "DCComponent.h"
#include "DCData.h"
#include <vector>
#include <utility>

/*CBMC graph of a branched/cyclic molecule
* The Decoupled/Coupled CBMC algorithm is represented by
* traversing a spanning tree of the graph.
*/

namespace cbmc
{

   class DCComponent;

   class DCGraph : public CBMC
   {
   public:
      DCGraph(System& sys, const Forcefield& ff,
         const MoleculeKind& kind, const Setup& set);

      void Build(TrialMol& oldMol, TrialMol& newMol, uint molIndex);
      ~DCGraph();

   private:
      struct Edge
      {
         uint destination;
         DCComponent* component;
         Edge(uint d, DCComponent* c) : destination(d), component(c) {}
      };
      struct Node
      {
         DCComponent* starting;
         std::vector<Edge> edges;
      };

      DCData data;
      std::vector<Node> nodes;
      std::vector<Edge> fringe;
      std::vector<bool> visited; //yes, I know
   };
}


#endif
