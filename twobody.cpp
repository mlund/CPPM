#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cylinder> Tspace;
typedef LennardJonesLB Tpairpot;

int main(int argc, char** argv) {
  InputMap mcp("twobody.json");

  Tspace spc(mcp);
  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::MassCenterConstrain<Tspace>(mcp, spc);

  Analysis::LineDistribution<> rdf(0.1);
  Analysis::CombinedAnalysis analysis(mcp, pot, spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  spc.load("state"); // load previous state, if any

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);
  while ( loop[0] ) {
    while ( loop[1] ) {
      mv.move();
      analysis.sample();
      rdf( spc.geo.dist( spc.groupList()[0]->cm, spc.groupList()[1]->cm ))++;

    } // end of micro loop

    cout << loop.timing();
    rdf.save("rdf.dat");

  } // end of macro loop

  cout << loop.info() + mv.info() + analysis.info() << endl;
}
