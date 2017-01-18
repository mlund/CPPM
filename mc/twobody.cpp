#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cylinder> Tspace;
//typedef Space<Geometry::Cuboid> Tspace;
typedef CombinedPairPotential<WeeksChandlerAndersen, Coulomb> Tpairpot;
//typedef CombinedPairPotential<HardSphere, Coulomb> Tpairpot;

int main(int argc, char** argv) {
  InputMap mcp("twobody.json");

  Tspace spc(mcp);
  spc.load("state"); // load previous state, if any

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp) + Energy::MassCenterConstrain<Tspace>(mcp, spc);

  //pot.first.pairpot.first.customParameters( mcp["customlj"] );

  Analysis::LineDistribution<> rdf(0.1);
  Analysis::CombinedAnalysis analysis(mcp, pot, spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  //analysis.add( MeanForce<Tspace>(
  //      { {"groups", {0,1} }, {"nstep", 10} },
  //      pot, spc) );

  cout << atom.info() << spc.info() << pot.info() << textio::header("MC Simulation Begins!");

  MCLoop loop(mcp);
  while ( loop[0] ) {
    while ( loop[1] ) {
      mv.move();
      //spc.groupList()[0]->setMassCenter(spc);
      //spc.groupList()[1]->setMassCenter(spc);
      analysis.sample();
      rdf( spc.geo.dist( spc.groupList()[0]->cm, spc.groupList()[1]->cm ))++;

    } // end of micro loop

    cout << loop.timing() << std::flush;
    rdf.save("rdf.dat");

  } // end of macro loop

  cout << loop.info() + mv.info() + analysis.info() << endl;
}
