#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

#ifdef CUBOID
typedef Space<Geometry::Cuboid> Tspace;
#else
typedef Space<Geometry::Cylinder> Tspace;
#endif
typedef CombinedPairPotential<WeeksChandlerAndersen, Coulomb> Tpairpot;

int main(int argc, char** argv) {
  InputMap mcp("twobody.json");

  Tspace spc(mcp);
  spc.load("state");

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp) + Energy::MassCenterConstrain<Tspace>(mcp, spc);

  Analysis::CombinedAnalysis analysis(mcp, pot, spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);
  while ( loop[0] ) {
    while ( loop[1] ) {
      mv.move();
      analysis.sample();
    } // end of micro loop

    cout << loop.timing() << std::flush;

  } // end of macro loop

  cout << loop.info() << mv.info() << analysis.info() << endl;
}
