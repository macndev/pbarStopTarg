// includes
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <algorithm>

#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGeneral.h"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "SeedService/inc/SeedService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/GenParticle.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/EventWeight.hh"
#include "Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Mu2eUtilities/inc/MuonCaptureSpectrum.hh"
#include "Mu2eUtilities/inc/SimpleSpectrum.hh"
#include "Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Mu2eUtilities/inc/Table.hh"
#include "Mu2eUtilities/inc/RootTreeSampler.hh"
#include "GeneralUtilities/inc/RSNTIO.hh"
#include "StoppingTargetGeom/inc/StoppingTarget.hh"
#include "StoppingTargetGeom/inc/TargetFoil.hh"

#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

namespace mu2e{

  class pbarStopTarg : public art::EDProducer {
    fhicl::ParameterSet psphys_;

    // double rhoInternal_;
    double elow_;
    double ehi_;

    BinnedSpectrum spectrum_;

    int verbosityLevel_;

    art::RandomNumberGenerator::base_engine_t& eng_;
    const double czmax_;
    const double czmin_;
    CLHEP::RandGeneral* randSpectrum_;
    RandomUnitSphere randUnitSphere_;
    RandomUnitSphere randUnitSphereExt_;
    CLHEP::RandFlat randFlat_;

    MuonCaptureSpectrum muonCaptureSpectrum_;

    bool doHistograms_;

    TH1F* _hmomentum;
    TH1F* _hCosz;
    TH1F* _hWeight;
    TH1F* _htZero;
    TH1F* _hxPos;
    TH1F* _hyPos;
    TH1F* _hzPos;
    
    // functions

    public:
      explicit pbarStopTarg(const fhicl::ParameterSet& pset);
      ~pbarStopTarg();
      virtual void produce(art::Event& event);
    };

    pbarStopTarg::pbarStopTarg(const fhicl::ParameterSet& pset)
      : EDProducer{pset}
      , psphys_             (pset.get<fhicl::ParameterSet>("physics"))
      // , rhoInternal_        (psphys_.get<double>("rhoInternal"))
      , spectrum_           (BinnedSpectrum(psphys_))
      , verbosityLevel_     (pset.get<int>("verbosityLevel", 0))
      , eng_                (createEngine(art::ServiceHandle<SeedService>()->getSeed()))
      , czmax_              (pset.get<double>("czmax",  1.))
      , czmin_              (pset.get<double>("czmin", -1.))
      , randUnitSphere_     (eng_)
      , randUnitSphereExt_  (eng_, czmax_, czmin_)
      , randFlat_           (eng_)
      , muonCaptureSpectrum_(&randFlat_,&randUnitSphere_)
	//   , stops_              (eng_, pset.get<fhicl::ParameterSet>("muonStops"))
      , doHistograms_       (pset.get<bool>("doHistograms",true ) )
    {
      produces<mu2e::GenParticleCollection>();
      produces<mu2e::EventWeight>();

      // fractionSpectrum_ = 0.;

      if(verbosityLevel_ > 0) {                              
        std::cout<<"pbarStopTarg: "
                 <<" stopped particles"
                 <<std::endl;
        std::cout<<"pbarStopTarg: producing pbar " << std::endl;
      }

      // double me_  = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::e_minus ).ref().mass().value();
      // double mmu_ = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::mu_minus).ref().mass().value();
      
      // double mpbar_ = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::anti_proton).ref().mass().value();

    randSpectrum_ = new CLHEP::RandGeneral(eng_, spectrum_.getPDF(), spectrum_.getNbins());

    if ( doHistograms_ ) {
      art::ServiceHandle<art::TFileService> tfs;
      art::TFileDirectory tfdir = tfs->mkdir( "pbarStopTarg" );

      _hmomentum       = tfdir.make<TH1F>("hmomentum", "Produced pbar momentum, RMC", 70,  0.,  140.  );
      _hCosz           = tfdir.make<TH1F>("hCosz", "Produced pbar cosz, RMC", 200,  -1.,  1.  );
      _htZero          = tfdir.make<TH1F>("htZero"         , "Stopped pbar time", 100,0.,2000.);
      _hWeight         = tfdir.make<TH1F>("hWeight"        , "Event Weight ", 100,0.,2.e-5);
      _hxPos           = tfdir.make<TH1F>("hxPos", "Stopped pbar x position", 100, 50., -50.);
      _hyPos           = tfdir.make<TH1F>("hyPos", "Stopped pbar y position", 100, 50., -50.);
      _hzPos           = tfdir.make<TH1F>("hzPos", "Stopper pbar z position",100, 0., 1000.);
    }
  }

    pbarStopTarg::~pbarStopTarg() {
      delete randSpectrum_;
   }

  // make structure to hold random x, y, z position in stopping target and time t
  struct stop {
    double x;
    double y;
    double z;
    double t;
  };

  stop eventVec;

  void pbarStopTarg::produce(art::Event& event) {

    std::unique_ptr<GenParticleCollection> output(new GenParticleCollection);

    // stopping target targ object, added stopping target include for this, get geom variables
    auto const& target = *GeomHandle<StoppingTarget>();

    unsigned int nfoils = target.nFoils();

    for(unsigned int i=0; i<nfoils; ++i)
      {
	const mu2e::TargetFoil& foilTarg = target.foil(i);
	int id = foilTarg.id();
	double cx = foilTarg.centerInDetectorSystem().x();
	double cy = foilTarg.centerInDetectorSystem().y();
	double cz = foilTarg.centerInDetectorSystem().z();
	double oradius = foilTarg.rOut();
	double halfThickness = foilTarg.halfThickness();
	double iradius = foilTarg.rIn();

	std::cout<<"The coordinates for foil "<<i<<" are: x center - "<<cx<<" y center - "<<cy<<" z center - "<<cz<<" inner radius - "<<iradius<<" outer radius - "<<oradius<<" and half thickness - "<<halfThickness<<" id - "<<id<<". Onto the next foil... "<<std::endl;
	  
      }
    
    // choose random foil and random position in x, y, z (mm)
    int n = randFlat_.fireInt(1, nfoils);
    
    //  double randrad = randFlat_.fire(foilTarg.iradius(n), foilTarg.oradius(n));
    double randrad = randFlat_.fire(target.foil(n).rIn(), target.foil(n).rOut());
    double randrho = randFlat_.fire(0, 2 * CLHEP::pi);

    eventVec.x = randrad * cos(randrho);
    eventVec.y = randrad * sin(randrho);
    eventVec.z = randFlat_.fire(target.foil(n).centerInDetectorSystem().z()-target.foil(n).halfThickness(), 
				target.foil(n).centerInDetectorSystem().z()+target.foil(n).halfThickness());

    // set time of ejection from stopping target, tZero; set x, y, z below
    eventVec.t = 1000.; // constant for now
    
    // define pos vector using random positions
    const CLHEP::Hep3Vector pos(eventVec.x, eventVec.y, eventVec.z);
    
    if (doHistograms_){
        _htZero->Fill(eventVec.t);
	_hxPos->Fill(eventVec.x);
	_hyPos->Fill(eventVec.y);
        _hzPos->Fill(eventVec.z);
    }

    double mpbar_ = GlobalConstantsHandle<ParticleDataTable>()->particle(PDGCode::anti_proton).ref().mass().value();

    double kLow = 49.9;
    double kHi = 50.0;
    double kinen = randFlat_.fire(kLow, kHi); // choose kinetic energy; matches elow, ehi
    
    double en = mpbar_*(CLHEP::c_light)*(CLHEP::c_light)*1000;  // rest energy (with conversion to MeV from GeV mass)
    
    double energy = kinen + en;  // kinetic + rest energy

    double mom = std::sqrt((energy * energy) - (en * en)); //momentum from energy

    double weight = 1.;  // weight = 1

    CLHEP::HepLorentzVector pbar(randUnitSphereExt_.fire(mom),energy);
    output->emplace_back( PDGCode::anti_proton,
                          GenId::pbarStopTarg,
                          pos,
                          pbar,
                          eventVec.t );

    event.put(std::move(output));

    //event.put(std::move(pw));

    if ( doHistograms_ ) {
      _hmomentum->Fill(mom);
      _hCosz->Fill(pbar.cosTheta());
      _hWeight->Fill(weight);
    }
   
    if (verbosityLevel_ > 0) {
      std::cout << "original pbar energy = " << energy << " and pbar mass = " << mpbar_ <<  std::endl;
      std::cout << "stop time = " << eventVec.t << std::endl;
    }
  }
}



 // namespace mu2e

DEFINE_ART_MODULE(mu2e::pbarStopTarg);
