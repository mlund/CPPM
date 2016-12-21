#include <faunus/faunus.h>

using namespace Faunus;
using namespace Faunus::Potential;

typedef Space<Geometry::Cylinder> Tspace;
typedef CombinedPairPotential<LennardJonesLB, Coulomb> Tpairpot;

/**
 * @brief Mean force between two groups
 *
 * This will average the mean force along the mass center
 * connection line between two groups. Currently, results
 * are output only via the `json()` function.
 *
 * JSON keys | Description
 * --------- | ---------------
 *  `groups` | Exactly two group index (array)
 *  `nsteps` | Sample interval
 */
template<class Tspace>
class MeanForce : public Analysis::AnalysisBase {

  typedef Energy::Energybase<Tspace> Tenergy;
  Tenergy& pot;
  Tspace& spc;
  Average<double> mf1, mf2;
  size_t g1, g2;

  void _sample() override {

    Point f1 = {0,0,0}; // net force on group 1
    Point f2 = {0,0,0}; // net force on group 2
    auto &g = spc.groupList(); // list of groups
    auto &p = spc.p;           // particle vector

    // force all others, k <-> g1 and g2
    for (size_t k=0; k!=g.size(); k++)
      if (k!=g1)
        if (k!=g2)
          for (auto i : *g[k]) {
            for (auto j : *g[g1])
              f1 += pot.f_p2p( p[i], p[j] );
            for (auto j : *g[g2])
              f2 += pot.f_p2p( p[i], p[j] );
          }
    // force g1<->g2
    for (auto i : *g[g2])
      for (auto j : *g[g1]) {
        Point f = pot.f_p2p( p[i], p[j] );
        f1 += f;
        f2 -= f;
      }

    // COM-COM unit vector and mean force
    Point r = spc.geo.vdist( g[g1]->cm, g[g2]->cm );
    mf1 += f1.dot( r/r.norm() );
    mf2 += f2.dot( -r/r.norm() );
  }

  string _info() override { return string(); }

  Tmjson _json() override
  {
    assert( mf1.cnt>0 && mf2.cnt>0 );
    return {
      { name,
        {
          { "groups", {g1, g2} },
          { "meanforce", { mf1.avg(), mf2.avg() } },
          { "forceunit", "kT/angstrom" }
        }
      }
    };
  }

  public:
  MeanForce( Tmjson j, Tenergy &pot, Tspace &spc ) : AnalysisBase(j), pot(pot), spc(spc)
  {
    name = "Mean force";
    g1 = g2 = -1;
    if (j["groups"].is_array()) {
      vector<size_t> v = j["groups"];
      if (v.size()==2) {
        g1 = v[0];
        g2 = v[1];
        if (g1>=0)
          if (g2>=0)
            if (g1!=g2)
              if (g1 < spc.groupList().size() )
                if (g2 < spc.groupList().size() )
                  return;
      }
    }
    throw std::runtime_error(name + ": group array must contain \
        exactly two distinct and valid group index");
  }

};

int main(int argc, char** argv) {
  InputMap mcp("twobody.json");

  Tspace spc(mcp);
  spc.load("state"); // load previous state, if any

  auto pot = Energy::Nonbonded<Tspace,Tpairpot>(mcp)
    + Energy::MassCenterConstrain<Tspace>(mcp, spc);

  pot.first.pairpot.first.customParameters( mcp["customlj"] );

  Analysis::LineDistribution<> rdf(0.1);
  Analysis::CombinedAnalysis analysis(mcp, pot, spc);
  Move::Propagator<Tspace> mv(mcp,pot,spc);

  analysis.add( MeanForce<Tspace>(
        { {"groups", {0,1} }, {"nstep", 20} },
        pot, spc) );

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
