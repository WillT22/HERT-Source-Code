///
/// Step: 1) Define the detectors (box --> logical volume (material) --> physical volume (locations etc))
///       2) SetSensDet (based on the material, it will identify which detector it is and set the sensitivity accordingly)
///       3) SetSensitiveDetector
/// Update: September 13, 2018
//          - Update the geometry for Reptile2 (include an outer diameters)
//  Revised by H. Zhao: Aug 12, 2020
//          - Change Be window thickness to 1.5 mm (median energy deposited on detector 1 inner of 1 MeV ~ 500 keV)
//          - Change detector thickness to 1.5/3 mm to match REPT (REPT R1 and R2 are 1.5 mm thick; R3 - R9 are comprised of two 1.5 mm thick silicon detectors)
//          - Reduce the distances between the detectors to reduce the influence from electron scattering
//  Revised by H. Zhao: Aug 14, 2020
//          - Change Be window thickness to 1.0 mm to make sure ~50% of 1 MeV electrons can reach R2 (and change the locations of tooth6 and Be window back)
//  HERT Design #4, revised by H. Zhao: Aug 26, 2020
//  - Be window: 1.5 mm
//  - Similar detector design R1 - R5 as REPT
//  - Longer collimater before Be window to achieve a smaller FOV (30 deg); remove embeded materials
//  - Increased shielding at the end: but still for boresight shooting only
//  - Add one more Si detector (as a whole; no inner/outer part) at the end of the stack of design #1

// ADD CHANGED NOTES

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "SetSensDet.hh"
#include <G4SubtractionSolid.hh>
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <G4NistManager.hh>
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"

DetectorConstruction::DetectorConstruction()
	:  //fDetectorMessenger(0),
	
	logic_w(0), logic_S1(0), logic_frontcoll(0), logic_coll(0), 
	logic_coll_embed1(0), logic_coll_embed2(0), logic_coll_embed3(0), logic_coll_embed4(0), logic_coll_embed5(0), 
	logic_coll_embed_Al1(0), logic_coll_embed_Al2(0), logic_coll_embed_Al3(0), logic_coll_embed_Al4(0), 
	logic_coll_tooth_W1(0), logic_coll_tooth_W2(0), logic_coll_tooth_W3(0), logic_coll_tooth_W4(0), logic_coll_tooth_W5(0), logic_coll_tooth_W6(0), logic_coll_tooth_W7(0), 
	logic_al_ann1(0), logic_al_chm(0), logic_al_end(0), logic_al_endann(0), logic_BeWin(0), 
	logic_w_ann(0), logic_w_ann1(0), logic_w_ann2(0), logic_w_ann3(0), logic_w_ann4(0), logic_w_ann_Al1(0), logic_w_ann_Al2(0), logic_w_ann_Al3(0), 
	logic_w_chm(0), logic_w_chm1(0), logic_w_chm2(0), logic_w_chm3(0), logic_w_chm4(0), logic_w_chm5(0), 
	logic_w_chm_Al1(0), logic_w_chm_Al2(0), logic_w_chm_Al3(0), logic_w_chm_Al4(0), 
	logic_w_end(0), logic_w_end1(0), logic_w_end2(0), logic_w_end3(0), logic_w_end4(0), logic_w_end5(0), 
	logic_w_end_Al1(0), logic_w_end_Al2(0), logic_w_end_Al3(0), logic_w_end_Al4(0), 
	logic_w_endenh(0), logic_w_endenh1(0), logic_w_endenh2(0), logic_w_endenh3(0), logic_w_endenh_Al1(0), logic_w_endenh_Al2(0), 
	logic_d1(0), logic_d2(0), logic_d3(0), logic_d4(0), logic_d5(0), logic_d6(0), logic_d7(0), logic_d8(0), logic_d9(0), 
	logic_fastener1(0), logic_fastener2(0), logic_fastener3(0), 
	physi_w(0), physi_S1(0), physi_frontcoll(0), physi_coll(0), 
	physi_coll_embed1(0), physi_coll_embed2(0), physi_coll_embed3(0), physi_coll_embed4(0),	physi_coll_embed5(0), 
	physi_coll_embed_Al1(0), physi_coll_embed_Al2(0), physi_coll_embed_Al3(0), physi_coll_embed_Al4(0), 
	physi_coll_tooth_W1(0), physi_coll_tooth_W2(0), physi_coll_tooth_W3(0), physi_coll_tooth_W4(0),	physi_coll_tooth_W5(0), physi_coll_tooth_W6(0), physi_coll_tooth_W7(0),	
	physi_al_ann1(0), physi_al_chm(0), physi_al_end(0), physi_al_endann(0), physi_be(0), 
	physi_w_ann(0), physi_w_ann1(0), physi_w_ann2(0), physi_w_ann3(0), physi_w_ann4(0), physi_w_ann_Al1(0), physi_w_ann_Al2(0), physi_w_ann_Al3(0), 
	physi_w_chm(0), physi_w_chm1(0), physi_w_chm2(0), physi_w_chm3(0), physi_w_chm4(0), physi_w_chm5(0), 
	physi_w_chm_Al1(0), physi_w_chm_Al2(0), physi_w_chm_Al3(0), physi_w_chm_Al4(0), 
	physi_w_end(0), physi_w_end1(0), physi_w_end2(0), physi_w_end3(0), physi_w_end4(0), physi_w_end5(0),
	physi_w_end_Al1(0), physi_w_end_Al2(0), physi_w_end_Al3(0), physi_w_end_Al4(0), 
	physi_w_endenh(0), physi_w_endenh1(0), physi_w_endenh2(0), physi_w_endenh3(0), physi_w_endenh_Al1(0), physi_w_endenh_Al2(0), 
	physi_d1(0), physi_d2(0), physi_d3(0), physi_d4(0), physi_d5(0), physi_d6(0), physi_d7(0), physi_d8(0), physi_d9(0), 
	physi_fastener1(0), physi_fastener2(0), physi_fastener3(0)
{

    //G4cout << "----> DetectorMessenger-Before." << G4endl;          
      fDetectorMessenger = new DetectorMessenger(this);
      windowDepth    = 1.5 *mm; //set window depth; to incorporate a thicker Be window, z of Be window has been adjusted (assuming the tungsten annulus is holding Be window)
      }
DetectorConstruction::~DetectorConstruction()
{
    //;
  delete fDetectorMessenger;
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
        //G4cout << "----> DetectorMessenger-After." << G4endl;  
	//____________________ materials ____________________ 

	G4double a;  // atomic mass
	G4double z;  // atomic number
	G4double density;
	G4int ncomponents;
	G4double fractionmass;
	
	// Use NIST database for elements and materials wherever possible.
	G4NistManager* man = G4NistManager::Instance();
	man->SetVerbose(1);
	G4Element* C  = man->FindOrBuildElement("C");
	G4Element* Mn = man->FindOrBuildElement("Mn");
	G4Element* Mg = man->FindOrBuildElement("Mg");
	G4Element* P = man->FindOrBuildElement("P");
	G4Element* S = man->FindOrBuildElement("S");
	G4Element* SiEle = man->FindOrBuildElement("Si");
	G4Element* Cr = man->FindOrBuildElement("Cr");
	G4Element* Fe = man->FindOrBuildElement("Fe");
	G4Element* Ni = man->FindOrBuildElement("Ni");
	G4Element* N = man->FindOrBuildElement("N");
	G4Element* Ti = man->FindOrBuildElement("Ti");
	G4Element* Zn = man->FindOrBuildElement("Zn");
	G4Element* Zr = man->FindOrBuildElement("Zr");
	G4Element* AlEle = man->FindOrBuildElement("Al");
	G4Element* elH = man->FindOrBuildElement("H");
	G4Element* elO = man->FindOrBuildElement("O");
	G4Element* elN = man->FindOrBuildElement("N");


	// Aluminum
	G4Material* Al = 
	new G4Material("Aluminum", z= 13.,  a= 26.9815*g/mole, density= 2.6989*g/cm3);

	// Tungsten
    G4Material* W =
        // spring '08 value: new G4Material("Tungsten", z= 74., a= 183.84*g/mole, density= 19.25*g/cm3); // wikipedia
	new G4Material("Tungsten", z= 69.5, a= 171.8*g/mole, density= 16.7*g/cm3); // 90% W; 10% Cu

    // Tantalum
	G4Material* Ta = new G4Material("Tantalum",z= 73., a= 180.95*g/mole, density= 16.65*g/cm3);

	// Silicon
	G4Material* Si = new G4Material("Silicon",z= 14., a= 28.0855*g/mole, density= 2.33*g/cm3);
  
	// Beryllium
	G4Material* Be = new G4Material("Beryllium",z= 4., a= 9.0122*g/mole, density= 1.848*g/cm3);
	
	// Aluminum (6061-T6) Reference: https://www.mcmaster.com/89015k111
	G4Material* Alalloy= new G4Material("AluminiumAlloy", density= 2.768*g/cm3, ncomponents=11);
	Alalloy->AddElement(AlEle, fractionmass=0.9525);
	Alalloy->AddElement(Cr, fractionmass=0.008);
	Alalloy->AddElement(C, fractionmass=0.004);
	Alalloy->AddElement(Fe, fractionmass=0.007);
    Alalloy->AddElement(Mg, fractionmass=0.012);
	Alalloy->AddElement(Mn, fractionmass=0.0015);
	Alalloy->AddElement(Ni, fractionmass=0.0005);
	Alalloy->AddElement(SiEle, fractionmass=0.008);
    Alalloy->AddElement(Ti, fractionmass=0.0015);
	Alalloy->AddElement(Zn, fractionmass=0.0025);
	Alalloy->AddElement(Zr, fractionmass=0.0025);
	

	// Stainless Steel - 304 (Based on http://www.worldstainless.org/Files/issf/non-image-files/PDF/Atlas_Grade_datasheet_-_all_datasheets_rev_Aug_2013.pdf) 
	G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=9);
	StainlessSteel->AddElement(C, fractionmass=0.0007);
	StainlessSteel->AddElement(SiEle, fractionmass=0.0075);
	StainlessSteel->AddElement(P, fractionmass=0.00045);
	StainlessSteel->AddElement(S, fractionmass=0.0003);
	StainlessSteel->AddElement(Cr, fractionmass=0.195);
	StainlessSteel->AddElement(Mn, fractionmass=0.02);
	StainlessSteel->AddElement(Fe, fractionmass=0.67);
	StainlessSteel->AddElement(Ni, fractionmass=0.105);
	StainlessSteel->AddElement(N, fractionmass=0.001);

	// Epoxy Base
	G4Material* EpoxyBase = new G4Material("EpoxyBase", density = 1.17 * g / cm3, ncomponents = 3);
	EpoxyBase->AddElement(elH, 21); // Hydrogen
	EpoxyBase->AddElement(elO, 3); // Oxygen
	EpoxyBase->AddElement(C, 28); // Carbon

	//Epoxy Accelerator
	G4Material* EpoxyAccelerator = new G4Material("EpoxyAccelerator", density = 0.97 * g / cm3, ncomponents = 4);
	EpoxyAccelerator->AddElement(elH, 24); // Hydrogen
	EpoxyAccelerator->AddElement(elN, 2); // Nitrogen
	EpoxyAccelerator->AddElement(elO, 3); // Oxygen
	EpoxyAccelerator->AddElement(C, 10); // Carbon
	//Density reference: https://www.globalspec.com/industrial-directory/density_epoxy_adhesives
	G4Material* Epoxy = new G4Material("Epoxy", density = 1.31 * g / cm3, ncomponents = 2);
	Epoxy->AddMaterial(EpoxyBase, 50 * perCent);
	Epoxy->AddMaterial(EpoxyAccelerator, 50 * perCent);

	//Density reference: https://www.sciencedirect.com/science/article/pii/S0168583X15004437
	G4Material* EpoxyTungsten = new G4Material("EpoxyW", density = 11 * g / cm3, ncomponents = 2);
	EpoxyTungsten->AddMaterial(Epoxy, 36 * perCent);
	EpoxyTungsten->AddMaterial(W, 64 * perCent);
	
	// Near Vacuum
	G4double atomicNumber = 1.;
	G4double massOfMole = 1.008*g/mole;
	G4double vac_density = 1.e-25*g/cm3;
	G4double temperature = 2.73*kelvin;
	G4double pressure = 3.e-18*pascal;
	G4Material* Vacuum = new G4Material("interGalactic", atomicNumber,massOfMole, vac_density, kStateGas,temperature, pressure);

	/*// air Reference: https://apc.u-paris.fr/~franco/g4doxy4.10/html/class_materials.html
	density = 1.2929e-03 * g / cm3;  // at 20 degree
	G4Material * Air = new G4Material("Air", density, nel = 3,
									   kStateGas, expTemp);
	G4double ttt = 75.47 + 23.20 + 1.28;
	Air->AddElement(elN, massfraction = 75.47 / ttt);
	Air->AddElement(elO, massfraction = 23.20 / ttt);
	Air->AddElement(elAr, massfraction = 1.28 / ttt);*/

	/*//Air?
	G4double atomicNumber = 7.;
	G4double massOfMole = 14*g/mole;
	G4double vac_density = 1.5e-3*g/cm3;
	G4double temperature = 296.*kelvin;
	G4double pressure = 101325.*pascal;
	G4Material* Vacuum = new G4Material("interGalactic", atomicNumber,massOfMole, vac_density, kStateGas,temperature, pressure);*/
		
	//Set the window material and depth
	windowMaterialBe = Be; // this is to ensure the Be and Al properties are the same as what we set and also limit the options we have for window material
	windowMaterialAl = Al;
	windowMaterial = windowMaterialBe;

	//____________________ solids ____________________

	// world volume (w)
	G4double w_fl = 40.0*cm; // full length
	G4double w_hl = 0.5*w_fl; // half length
	
	G4Box* solid_w = new G4Box("world",w_hl,w_hl,w_hl);
	logic_w = new G4LogicalVolume(solid_w,Vacuum,"world",0,0,0);
	physi_w = new G4PVPlacement(0,G4ThreeVector(),logic_w,"world",0,false,0);

	G4double gapR = 1.08 * mm; // gap from back of detect to front of next detector
	G4double gapB = 1.0 * mm; // gap from back of detect to front of next detector

	//Source Spherical Cap
	G4double pRmin = 85;// *mm;
	G4double pRmax = pRmin + 0.2 * mm;
	G4double pSPhi = 0 * deg;
	G4double pDPhi = 360 * deg;
	G4double pSTheta = 165 * deg; //165 for 15 degree cap
	G4double pDTheta = 180 * deg - pSTheta;
	G4double S1_x = 0.0 * mm; // x location
	G4double S1_y = 0.0 * mm; // y location
	G4double S1_z = 0.5 * 1.5 * mm; // z location-centered on first detector  d1_z
	const G4String& pName = "source";

	G4Sphere* S1 = new G4Sphere(pName, pRmin, pRmax, pSPhi, pDPhi, pSTheta, pDTheta);
	logic_S1 = new G4LogicalVolume(S1, Vacuum, "source", 0, 0, 0);
	physi_S1 = new G4PVPlacement(0, G4ThreeVector(S1_x, S1_y, S1_z), logic_S1, "source", logic_w, false, 0);

	// detectors (d)
	// detector 1 (d1)
	G4double d1_d = 1.5 * mm; // depth
	G4double d1_hd = 0.5 * d1_d * mm; // half depth
	G4double d1_ir = 0.0 * mm; // inner radius
	G4double d1_or = 10.0 * mm; // outer radius
	G4double d1_sta = 0.0 * deg; // start angle
	G4double d1_spa = 360 * deg; // span angle
	G4double d1_x = 0.0 * mm; // x location
	G4double d1_y = 0.0 * mm; // y location
	G4double d1_z = d1_hd; // z location
	G4Tubs* solid_d1 = new G4Tubs("detector_1", d1_ir, d1_or, d1_hd, d1_sta, d1_spa);
	logic_d1 = new G4LogicalVolume(solid_d1, Si, "detector_1", 0, 0, 0);
	physi_d1 = new G4PVPlacement(0, G4ThreeVector(d1_x, d1_y, d1_z), logic_d1, "detector_1", logic_w, false, 0);

	// detector 2 (d2)
	G4double d2_d = 1.5 * mm; // depth
	G4double d2_hd = 0.5 * d2_d * mm; // half depth
	G4double d2_ir = 0.0 * mm; // inner radius
	G4double d2_or = 20.0 * mm; // outer radius
	G4double d2_sta = 0.0 * deg; // start angle
	G4double d2_spa = 360 * deg; // span angle
	G4double d2_x = 0.0 * mm; // x location
	G4double d2_y = 0.0 * mm; // y location
	G4double d2_z = d1_z + d2_d + gapR; // z location

	G4Tubs* solid_d2 = new G4Tubs("detector_2", d2_ir, d2_or, d2_hd, d2_sta, d2_spa);
	logic_d2 = new G4LogicalVolume(solid_d2, Si, "detector_2", 0, 0, 0);
	physi_d2 = new G4PVPlacement(0, G4ThreeVector(d2_x, d2_y, d2_z), logic_d2, "detector_2", logic_w, false, 0);

	// detector 3 (d3)
	G4double d3_d = 1.5 * mm; // depth
	G4double d3_hd = 0.5 * d3_d * mm; // half depth
	G4double d3_ir = 0.0 * mm; // inner radius
	G4double d3_or = 20.0 * mm; // outer radius
	G4double d3_sta = 0.0 * deg; // start angle
	G4double d3_spa = 360 * deg; // span angle
	G4double d3_x = 0.0 * mm; // x location
	G4double d3_y = 0.0 * mm; // y location
	G4double d3_z = d2_z + d3_d + gapB; // z location

	G4Tubs* solid_d3 = new G4Tubs("detector_3", d3_ir, d3_or, d3_hd, d3_sta, d3_spa);
	logic_d3 = new G4LogicalVolume(solid_d3, Si, "detector_3", 0, 0, 0);
	physi_d3 = new G4PVPlacement(0, G4ThreeVector(d3_x, d3_y, d3_z), logic_d3, "detector_3", logic_w, false, 0);

	// detector 4 (d4)
	G4double d4_d = 1.5 * mm; // depth
	G4double d4_hd = 0.5 * d4_d * mm; // half depth
	G4double d4_ir = 0.0 * mm; // inner radius
	G4double d4_or = 20.0 * mm; // outer radius
	G4double d4_sta = 0.0 * deg; // start angle
	G4double d4_spa = 360 * deg; // span angle
	G4double d4_x = 0.0 * mm; // x location
	G4double d4_y = 0.0 * mm; // y location
	G4double d4_z = d3_z + d4_d + gapR; // z location

	G4Tubs* solid_d4 = new G4Tubs("detector_4", d4_ir, d4_or, d4_hd, d4_sta, d4_spa);
	logic_d4 = new G4LogicalVolume(solid_d4, Si, "detector_4", 0, 0, 0);
	physi_d4 = new G4PVPlacement(0, G4ThreeVector(d4_x, d4_y, d4_z), logic_d4, "detector_4", logic_w, false, 0);

	// detector 5 (d5)
	G4double d5_d = 1.5 * mm; // depth
	G4double d5_hd = 0.5 * d5_d * mm; // half depth
	G4double d5_ir = 0.0 * mm; // inner radius
	G4double d5_or = 20.0 * mm; // outer radius
	G4double d5_sta = 0.0 * deg; // start angle
	G4double d5_spa = 360 * deg; // span angle
	G4double d5_x = 0.0 * mm; // x location
	G4double d5_y = 0.0 * mm; // y location
	G4double d5_z = d4_z + d5_d + gapB; // z location

	G4Tubs* solid_d5 = new G4Tubs("detector_5", d5_ir, d5_or, d5_hd, d5_sta, d5_spa);
	logic_d5 = new G4LogicalVolume(solid_d5, Si, "detector_5", 0, 0, 0);
	physi_d5 = new G4PVPlacement(0, G4ThreeVector(d5_x, d5_y, d5_z), logic_d5, "detector_5", logic_w, false, 0);

	// detector 6 (d6)
	G4double d6_d = 1.5 * mm; // depth
	G4double d6_hd = 0.5 * d6_d * mm; // half depth
	G4double d6_ir = 0.0 * mm; // inner radius
	G4double d6_or = 20.0 * mm; // outer radius
	G4double d6_sta = 0.0 * deg; // start angle
	G4double d6_spa = 360 * deg; // span angle
	G4double d6_x = 0.0 * mm; // x location
	G4double d6_y = 0.0 * mm; // y location
	G4double d6_z = d5_z + d6_d + gapR; // z location

	G4Tubs* solid_d6 = new G4Tubs("detector_6", d6_ir, d6_or, d6_hd, d6_sta, d6_spa);
	logic_d6 = new G4LogicalVolume(solid_d6, Si, "detector_6", 0, 0, 0);
	physi_d6 = new G4PVPlacement(0, G4ThreeVector(d6_x, d6_y, d6_z), logic_d6, "detector_6", logic_w, false, 0);

	// detector 7 (d7)
	G4double d7_d = 1.5 * mm; // depth
	G4double d7_hd = 0.5 * d7_d * mm; // half depth
	G4double d7_ir = 0.0 * mm; // inner radius
	G4double d7_or = 20.0 * mm; // outer radius
	G4double d7_sta = 0.0 * deg; // start angle
	G4double d7_spa = 360 * deg; // span angle
	G4double d7_x = 0.0 * mm; // x location
	G4double d7_y = 0.0 * mm; // y location
	G4double d7_z = d6_z + d7_d + gapB; // z location

	G4Tubs* solid_d7 = new G4Tubs("detector_7", d7_ir, d7_or, d7_hd, d7_sta, d7_spa);
	logic_d7 = new G4LogicalVolume(solid_d7, Si, "detector_7", 0, 0, 0);
	physi_d7 = new G4PVPlacement(0, G4ThreeVector(d7_x, d7_y, d7_z), logic_d7, "detector_7", logic_w, false, 0);

	// detector 8 (d8)
	G4double d8_d = 1.5 * mm; // depth
	G4double d8_hd = 0.5 * d8_d * mm; // half depth
	G4double d8_ir = 0.0 * mm; // inner radius
	G4double d8_or = 20.0 * mm; // outer radius
	G4double d8_sta = 0.0 * deg; // start angle
	G4double d8_spa = 360 * deg; // span angle
	G4double d8_x = 0.0 * mm; // x location
	G4double d8_y = 0.0 * mm; // y location
	G4double d8_z = d7_z + d8_d + gapR; // z location

	G4Tubs* solid_d8 = new G4Tubs("detector_8", d8_ir, d8_or, d8_hd, d8_sta, d8_spa);
	logic_d8 = new G4LogicalVolume(solid_d8, Si, "detector_8", 0, 0, 0);
	physi_d8 = new G4PVPlacement(0, G4ThreeVector(d8_x, d8_y, d8_z), logic_d8, "detector_8", logic_w, false, 0);

	// detector 9 (d9)
	G4double d9_d = 1.5 * mm; // depth
	G4double d9_hd = 0.5 * d9_d * mm; // half depth
	G4double d9_ir = 0.0 * mm; // inner radius
	G4double d9_or = 20.0 * mm; // outer radius
	G4double d9_sta = 0.0 * deg; // start angle
	G4double d9_spa = 360 * deg; // span angle
	G4double d9_x = 0.0 * mm; // x location
	G4double d9_y = 0.0 * mm; // y location
	G4double d9_z = d8_z + d9_d + gapB; // z location

	G4Tubs* solid_d9 = new G4Tubs("detector_9", d9_ir, d9_or, d9_hd, d9_sta, d9_spa);
	logic_d9 = new G4LogicalVolume(solid_d9, Si, "detector_9", 0, 0, 0);
	physi_d9 = new G4PVPlacement(0, G4ThreeVector(d9_x, d9_y, d9_z), logic_d9, "detector_9", logic_w, false, 0);
	
	G4RotationMatrix* rotm_AlHousing = new G4RotationMatrix();
	rotm_AlHousing->rotateX(-90. * deg);
	rotm_AlHousing->rotateY(-30. * deg);

	//Rotation Matrices
	G4RotationMatrix* rotm_negZ = new G4RotationMatrix();
	rotm_negZ->rotateZ(-90. * deg);

	G4RotationMatrix* rotm_negY = new G4RotationMatrix();
	rotm_negY->rotateY(-90. * deg);

	G4RotationMatrix* rotm_Y = new G4RotationMatrix();
	rotm_Y->rotateY(90. * deg);

	G4RotationMatrix* rotm_X = new G4RotationMatrix();
	rotm_X->rotateX(90. * deg);

	G4RotationMatrix* rotm_negX = new G4RotationMatrix();
	rotm_negX->rotateX(-90. * deg);

	// aluminum front collimator
	G4double frontcoll_d = 1.00 * mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double frontcoll_hd = 0.5 * frontcoll_d * mm; // half depth
	G4double frontcoll_ir = 9.00 * mm; // inner radius
	G4double frontcoll_or = 23.5 * mm; // outer radius | (5 mm thick Al collimator)
	G4double frontcoll_sta = 0.0 * deg; // start angle
	G4double frontcoll_spa = 360 * deg; // span angle
	G4double frontcoll_x = 0.0 * mm; // x location
	G4double frontcoll_y = 0.0 * mm; // y location
	G4double frontcoll_z = frontcoll_hd; // z location
	
	G4Tubs* solid_frontcoll = new G4Tubs("Al_frontcollimator", frontcoll_ir, frontcoll_or, frontcoll_hd, frontcoll_sta, frontcoll_spa);
	logic_frontcoll = new G4LogicalVolume(solid_frontcoll, Al, "Al_frontcollimator", 0, 0, 0);
	physi_frontcoll = new G4PVPlacement(0, G4ThreeVector(frontcoll_x, frontcoll_y, frontcoll_z), logic_frontcoll, "Al_frontcollimator", logic_w, false, 0);

	// aluminum collimator
	G4double coll_d = 47 * mm; // depth
	G4double coll_hd = 0.5 * coll_d * mm; // half depth
	G4double coll_ir = 17.125 * mm; // inner radius
	G4double coll_or = 23.5 * mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta = 0.0 * deg; // start angle
	G4double coll_spa = 360 * deg; // span angle
	G4double coll_x = 0.0 * mm; // x location
	G4double coll_y = 0.0 * mm; // y location
	G4double coll_z = frontcoll_d + coll_hd; // z location
	
	G4Tubs* solid_coll = new G4Tubs("Al_collimator", coll_ir, coll_or, coll_hd, coll_sta, coll_spa);
	logic_coll = new G4LogicalVolume(solid_coll, Al, "Al_collimator", 0, 0, 0);
	physi_coll = new G4PVPlacement(0, G4ThreeVector(coll_x, coll_y, coll_z), logic_coll, "Al_collimator", logic_w, false, 0);

	// heavy shielding embedded in collimator 
	G4double coll_d_embed5 = 58.5 * mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed5 = 0.5 * coll_d_embed5 * mm; // half depth
	G4double coll_ir_embed5 = 15.0 * mm; // inner radius
	G4double coll_or_embed5 = 17.0 * mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed5 = 0.0 * deg; // start angle
	G4double coll_spa_embed5 = 360 * deg; // span angle
	G4double coll_x_embed5 = 0.0 * mm; // x location
	G4double coll_y_embed5 = 0.0 * mm; // y location
	G4double coll_z_embed5 = frontcoll_d + coll_hd_embed5; // z location
	
	G4Tubs* solid_coll_embed5 = new G4Tubs("embed_collimator5", coll_ir_embed5, coll_or_embed5, coll_hd_embed5, coll_sta_embed5, coll_spa_embed5);
	logic_coll_embed5 = new G4LogicalVolume(solid_coll_embed5, Ta, "embed_collimator5", 0, 0, 0);
	physi_coll_embed5 = new G4PVPlacement(0, G4ThreeVector(coll_x_embed5, coll_y_embed5, coll_z_embed5), logic_coll_embed5, "embed_collimator5", logic_w, false, 0);

	// collimator tooth - W1
	G4double coll_d_tooth_W1 = 1.0 * mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W1 = 0.5 * coll_d_tooth_W1 * mm; // half depth
	G4double coll_ir_tooth_W1 = 9.00 * mm; // inner radius
	G4double coll_or_tooth_W1 = 15.00 * mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W1 = 0.0 * deg; // start angle
	G4double coll_spa_tooth_W1 = 360 * deg; // span angle
	G4double coll_x_tooth_W1 = 0.0 * mm; // x location
	G4double coll_y_tooth_W1 = 0.0 * mm; // y location
	G4double coll_z_tooth_W1 = frontcoll_d + coll_hd_tooth_W1; // z location
	
	G4Tubs* solid_coll_tooth_W1 = new G4Tubs("collimator_tooth_W1", coll_ir_tooth_W1, coll_or_tooth_W1,coll_hd_tooth_W1, coll_sta_tooth_W1, coll_spa_tooth_W1);
	logic_coll_tooth_W1 = new G4LogicalVolume(solid_coll_tooth_W1, Ta, "collimator_tooth_W1", 0, 0, 0);
	physi_coll_tooth_W1 = new G4PVPlacement(0, G4ThreeVector(coll_x_tooth_W1, coll_y_tooth_W1, coll_z_tooth_W1),logic_coll_tooth_W1, "collimator_tooth_W1", logic_w, false, 0);
	
	// collimator tooth - W2
	G4double coll_d_tooth_W2 = 1.0 * mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W2 = 0.5 * coll_d_tooth_W2 * mm; // half depth
	G4double coll_ir_tooth_W2 = 9.00 * mm; // inner radius
	G4double coll_or_tooth_W2 = 15.00 * mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W2 = 0.0 * deg; // start angle
	G4double coll_spa_tooth_W2 = 360 * deg; // span angle
	G4double coll_x_tooth_W2 = 0.0 * mm; // x location
	G4double coll_y_tooth_W2 = 0.0 * mm; // y location
	G4double coll_z_tooth_W2 = frontcoll_d + 15.0 * mm + coll_hd_tooth_W2; // z location
	
	G4Tubs* solid_coll_tooth_W2 = new G4Tubs("collimator_tooth_W2", coll_ir_tooth_W2, coll_or_tooth_W2, coll_hd_tooth_W2, coll_sta_tooth_W2, coll_spa_tooth_W2);
	logic_coll_tooth_W2 = new G4LogicalVolume(solid_coll_tooth_W2, Ta, "collimator_tooth_W2", 0, 0, 0);
	physi_coll_tooth_W2 = new G4PVPlacement(0, G4ThreeVector(coll_x_tooth_W2, coll_y_tooth_W2, coll_z_tooth_W2),logic_coll_tooth_W2, "collimator_tooth_W2", logic_w, false, 0);
	
	// collimator tooth - W3
	G4double coll_d_tooth_W3 = 1.0 * mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W3 = 0.5 * coll_d_tooth_W3 * mm; // half depth
	G4double coll_ir_tooth_W3 = 9.00 * mm; // inner radius
	G4double coll_or_tooth_W3 = 15.00 * mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W3 = 0.0 * deg; // start angle
	G4double coll_spa_tooth_W3 = 360 * deg; // span angle
	G4double coll_x_tooth_W3 = 0.0 * mm; // x location
	G4double coll_y_tooth_W3 = 0.0 * mm; // y location
	G4double coll_z_tooth_W3 = frontcoll_d + 30 * mm + coll_hd_tooth_W3; // z location
	
	G4Tubs* solid_coll_tooth_W3 = new G4Tubs("collimator_tooth_W3", coll_ir_tooth_W3, coll_or_tooth_W3, coll_hd_tooth_W3, coll_sta_tooth_W3, coll_spa_tooth_W3);
	logic_coll_tooth_W3 = new G4LogicalVolume(solid_coll_tooth_W3, Ta, "collimator_tooth_W3", 0, 0, 0);
	physi_coll_tooth_W3 = new G4PVPlacement(0, G4ThreeVector(coll_x_tooth_W3, coll_y_tooth_W3, coll_z_tooth_W3), logic_coll_tooth_W3, "collimator_tooth_W3", logic_w, false, 0);
	
	// collimator tooth - W4
	G4double coll_d_tooth_W4 = 1.0 * mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W4 = 0.5 * coll_d_tooth_W4 * mm; // half depth
	G4double coll_ir_tooth_W4 = 9.00 * mm; // inner radius
	G4double coll_or_tooth_W4 = 15.00 * mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W4 = 0.0 * deg; // start angle
	G4double coll_spa_tooth_W4 = 360 * deg; // span angle
	G4double coll_x_tooth_W4 = 0.0 * mm; // x location
	G4double coll_y_tooth_W4 = 0.0 * mm; // y location
	G4double coll_z_tooth_W4 = frontcoll_d + 45 * mm + coll_hd_tooth_W4; // z location
	
	G4Tubs* solid_coll_tooth_W4 = new G4Tubs("collimator_tooth_W4", coll_ir_tooth_W4, coll_or_tooth_W4, coll_hd_tooth_W4, coll_sta_tooth_W4, coll_spa_tooth_W4);
	logic_coll_tooth_W4 = new G4LogicalVolume(solid_coll_tooth_W4, Ta, "collimator_tooth_W4", 0, 0, 0);
	physi_coll_tooth_W4 = new G4PVPlacement(0, G4ThreeVector(coll_x_tooth_W4, coll_y_tooth_W4, coll_z_tooth_W4), logic_coll_tooth_W4, "collimator_tooth_W4", logic_w, false, 0);
	
	// collimator tooth - W5
	G4double coll_d_tooth_W5 = 1.0 * mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W5 = 0.5 * coll_d_tooth_W5 * mm; // half depth
	G4double coll_ir_tooth_W5 = 9.00 * mm; // inner radius
	G4double coll_or_tooth_W5 = 17.00 * mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W5 = 0.0 * deg; // start angle
	G4double coll_spa_tooth_W5 = 360 * deg; // span angle
	G4double coll_x_tooth_W5 = 0.0 * mm; // x location
	G4double coll_y_tooth_W5 = 0.0 * mm; // y location
	G4double coll_z_tooth_W5 = frontcoll_d + 60 * mm + coll_hd_tooth_W5; // z location
	
	G4Tubs* solid_coll_tooth_W5 = new G4Tubs("collimator_tooth_W5", coll_ir_tooth_W5, coll_or_tooth_W5, coll_hd_tooth_W5, coll_sta_tooth_W5, coll_spa_tooth_W5);
	logic_coll_tooth_W5 = new G4LogicalVolume(solid_coll_tooth_W5, Ta, "collimator_tooth_W5", 0, 0, 0);
	physi_coll_tooth_W5 = new G4PVPlacement(0, G4ThreeVector(coll_x_tooth_W5, coll_y_tooth_W5, coll_z_tooth_W5), logic_coll_tooth_W5, "collimator_tooth_W5", logic_w, false, 0);

	// aluminum front annulus 1 (larger of the two)
	G4double al_ann1_d = 10 * mm; // depth
	G4double al_ann1_hd = 0.5 * al_ann1_d * mm; // half depth
	G4double al_ann1_ir = 17.0 * mm; // inner radius
	G4double al_ann1_or = 40 * mm; // outer radius
	G4double al_ann1_sta = 0.0 * deg; // start angle
	G4double al_ann1_spa = 360 * deg; // span angle
	G4double al_ann1_x = 0.0 * mm; // x location
	G4double al_ann1_y = 0.0 * mm; // y location
	G4double al_ann1_z = frontcoll_d + coll_d + al_ann1_hd; // z location
	
	G4Tubs* solid_al_ann1 = new G4Tubs("Al_annulus_1", al_ann1_ir, al_ann1_or, al_ann1_hd, al_ann1_sta, al_ann1_spa);
	logic_al_ann1 = new G4LogicalVolume(solid_al_ann1, Al, "Al_annulus_1", 0, 0, 0);
	physi_al_ann1 = new G4PVPlacement(0, G4ThreeVector(al_ann1_x, al_ann1_y, al_ann1_z), logic_al_ann1, "Al_annulus_1", logic_w, false, 0);

	// aluminum chamber
	G4double al_chm_d = 36.5 * mm; // depth: revised for HERT1
	G4double al_chm_hd = 0.5 * al_chm_d * mm; // half depth
	G4double al_chm_ir = 35 * mm; // inner radius | 30.5625 - 5.0 - 3.5 = 22.0625
	G4double al_chm_or = 40 * mm; // outer radius
	G4double al_chm_sta = 0.0 * deg; // start angle
	G4double al_chm_spa = 360 * deg; // span angle
	G4double al_chm_x = 0.0 * mm; // x location
	G4double al_chm_y = 0.0 * mm; // y location
	G4double al_chm_z = frontcoll_d + coll_d + al_ann1_d + al_chm_hd; // z location. Add 0.75 for width of delrin lip
	
	G4Tubs* solid_al_chm = new G4Tubs("Al_chamber", al_chm_ir, al_chm_or, al_chm_hd, al_chm_sta, al_chm_spa);
	logic_al_chm = new G4LogicalVolume(solid_al_chm, Al, "Al_chamber", 0, 0, 0);
	physi_al_chm = new G4PVPlacement(0, G4ThreeVector(al_chm_x, al_chm_y, al_chm_z), logic_al_chm, "Al_chamber", logic_w, false, 0);

	// aluminum end plate
	G4double al_end_d = 5.0 * mm; // depth
	G4double al_end_hd = 0.5 * al_end_d * mm; // half depth
	G4double al_end_ir = 0.0 * mm; // inner radius
	G4double al_end_or = 40 * mm; // outer radius
	G4double al_end_sta = 0.0 * deg; // start angle
	G4double al_end_spa = 360 * deg; // span angle
	G4double al_end_x = 0.0 * mm; // x location
	G4double al_end_y = 0.0 * mm; // y location
	G4double al_end_z = frontcoll_d + coll_d + al_ann1_d + al_chm_d + al_end_hd; // z location
	
	G4Tubs* solid_al_end = new G4Tubs("Al_end_plate", al_end_ir, al_end_or, al_end_hd, al_end_sta, al_end_spa);
	logic_al_end = new G4LogicalVolume(solid_al_end, Al, "Al_end_plate", 0, 0, 0);
	physi_al_end = new G4PVPlacement(0, G4ThreeVector(al_end_x, al_end_y, al_end_z), logic_al_end, "Al_end_plate", logic_w, false, 0);

	// beryllium window
	G4double be_d = windowDepth; // depth
	G4double be_hd = 0.5 * windowDepth * mm; // half depth
	G4double be_ir = 0.0 * mm; // inner radius
	G4double be_or = 17.0 * mm; // outer radius
	G4double be_sta = 0.0 * deg; // start angle
	G4double be_spa = 360 * deg; // span angle
	G4double be_x = 0.0 * mm; // x location
	G4double be_y = 0.0 * mm; // y location
	G4double be_z = frontcoll_d + 58.5 * mm + 0.5 * windowDepth; // z location

	G4Tubs* solid_be = new G4Tubs("FOV_disc", be_ir, be_or, 0.5 * windowDepth, be_sta, be_spa);
	logic_BeWin = new G4LogicalVolume(solid_be, windowMaterial, "FOV_disc", 0, 0, 0);
	physi_be = new G4PVPlacement(0, G4ThreeVector(be_x, be_y, be_z), logic_BeWin, "FOV_disc", logic_w, false, 0);

	G4cout << "----> The window's z is " << be_z << G4endl;
	G4cout << "----> The window's material is " << logic_BeWin->GetMaterial() << G4endl;
	G4cout << "----> The window's depth is " << windowDepth << "mm" << G4endl;

	// tungsten front annulus part 1
	G4double w_ann_d = 3.5 * mm; // depth
	G4double w_ann_hd = 0.5 * w_ann_d * mm; // half depth
	G4double w_ann_ir = 17.125 * mm; // inner radius
	G4double w_ann_or = 34.75 * mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta = 0.0 * deg; // start angle
	G4double w_ann_spa = 360 * deg; // span angle
	G4double w_ann_x = 0.0 * mm; // x location
	G4double w_ann_y = 0.0 * mm; // y location
	G4double w_ann_z = frontcoll_d + coll_d + al_ann1_d + w_ann_hd; // z location
	
	G4Tubs* solid_w_ann = new G4Tubs("W_annulus", w_ann_ir, w_ann_or, w_ann_hd, w_ann_sta, w_ann_spa);
	logic_w_ann = new G4LogicalVolume(solid_w_ann, W, "W_annulus", 0, 0, 0);
	physi_w_ann = new G4PVPlacement(0, G4ThreeVector(w_ann_x, w_ann_y, w_ann_z), logic_w_ann, "W_annulus", logic_w, false, 0);

	// Al Insert: Between tungsten part 1 and part 2 front annulus part 1
	G4double w_ann_d_Al1 = 0.5 * mm; // depth
	G4double w_ann_hd_Al1 = 0.5 * w_ann_d_Al1 * mm; // half depth
	G4double w_ann_ir_Al1 = 17.125 * mm; // inner radius
	G4double w_ann_or_Al1 = 34.75 * mm; // outer radius | 
	G4double w_ann_sta_Al1 = 0.0 * deg; // start angle
	G4double w_ann_spa_Al1 = 360 * deg; // span angle
	G4double w_ann_x_Al1 = 0.0 * mm; // x location
	G4double w_ann_y_Al1 = 0.0 * mm; // y location
	G4double w_ann_z_Al1 = frontcoll_d + coll_d + al_ann1_d + w_ann_d + w_ann_hd_Al1; // z location
	
	G4Tubs* solid_w_ann_Al1 = new G4Tubs("W_annulus_Al1", w_ann_ir_Al1, w_ann_or_Al1, w_ann_hd_Al1, w_ann_sta_Al1, w_ann_spa_Al1);
	logic_w_ann_Al1 = new G4LogicalVolume(solid_w_ann_Al1, Alalloy, "W_annulus_Al1", 0, 0, 0);
	physi_w_ann_Al1 = new G4PVPlacement(0, G4ThreeVector(w_ann_x_Al1, w_ann_y_Al1, w_ann_z_Al1), logic_w_ann_Al1, "W_annulus_Al1", logic_w, false, 0);

	// tungsten front annulus part 2
	G4double w_ann_d1 = 1.5 * mm; // depth
	G4double w_ann_hd1 = 0.5 * w_ann_d1 * mm; // half depth
	G4double w_ann_ir1 = 11.5 * mm; // inner radius
	G4double w_ann_or1 = 34.75 * mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta1 = 0.0 * deg; // start angle
	G4double w_ann_spa1 = 360 * deg; // span angle
	G4double w_ann_x1 = 0.0 * mm; // x location
	G4double w_ann_y1 = 0.0 * mm; // y location
	G4double w_ann_z1 = frontcoll_d + coll_d + al_ann1_d + w_ann_d + w_ann_d_Al1 + w_ann_hd1; // z location
	
	G4Tubs* solid_w_ann1 = new G4Tubs("W_annulus1", w_ann_ir1, w_ann_or1, w_ann_hd1, w_ann_sta1, w_ann_spa1);
	logic_w_ann1 = new G4LogicalVolume(solid_w_ann1, W, "W_annulus1", 0, 0, 0);
	physi_w_ann1 = new G4PVPlacement(0, G4ThreeVector(w_ann_x1, w_ann_y1, w_ann_z1), logic_w_ann1, "W_annulus1", logic_w, false, 0);

	// tungsten chamber
	G4double w_chm_d = 24.86 * mm; // depth
	G4double w_chm_hd = 0.5 * w_chm_d * mm; // half depth
	G4double w_chm_ir = 30.75 * mm; // inner radius
	G4double w_chm_or = 34.75 * mm; // outer radius
	G4double w_chm_sta = 0.0 * deg; // start angle
	G4double w_chm_spa = 360 * deg; // span angle
	G4double w_chm_x = 0.0 * mm; // x location
	G4double w_chm_y = 0.0 * mm; // y location
	G4double w_chm_z = frontcoll_d + coll_d + al_ann1_d + w_ann_d + w_ann_d_Al1 + w_ann_d1 + w_chm_hd; // z location
	
	G4Tubs* solid_w_chm = new G4Tubs("W_chamber", w_chm_ir, w_chm_or, w_chm_hd, w_chm_sta, w_chm_spa);
	logic_w_chm = new G4LogicalVolume(solid_w_chm, EpoxyTungsten, "W_chamber", 0, 0, 0);
	physi_w_chm = new G4PVPlacement(0, G4ThreeVector(w_chm_x, w_chm_y, w_chm_z), logic_w_chm, "W_chamber", logic_w, false, 0);

	// tungsten end plate
	G4double w_end_d = 5.0 * mm; // depth
	G4double w_end_hd = 0.5 * w_end_d * mm; // half depth
	G4double w_end_ir = 0.0 * mm; // inner radius
	G4double w_end_or = 34.75 * mm; // outer radius
	G4double w_end_sta = 0.0 * deg; // start angle
	G4double w_end_spa = 360 * deg; // span angle
	G4double w_end_x = 0.0 * mm; // x location
	G4double w_end_y = 0.0 * mm; // y location
	G4double w_end_z = frontcoll_d + coll_d + al_ann1_d + w_ann_d + w_ann_d_Al1 + w_ann_d1 + w_chm_d + w_end_hd; // z location
	
	G4Tubs* solid_w_end = new G4Tubs("W_end_plate", w_end_ir, w_end_or, w_end_hd, w_end_sta, w_end_spa);
	logic_w_end = new G4LogicalVolume(solid_w_end, W, "W_end_plate", 0, 0, 0);
	physi_w_end = new G4PVPlacement(0, G4ThreeVector(w_end_x, w_end_y, w_end_z), logic_w_end, "W_end_plate", logic_w, false, 0);

	// Al Insert: Shims placed behind detector stack
	G4double w_end_d_Al1 = 2 * mm; // depth
	G4double w_end_hd_Al1 = 0.5 * w_end_d_Al1 * mm; // half depth
	G4double w_end_ir_Al1 = 0.0 * mm; // inner radius
	G4double w_end_or_Al1 = 29.3 * mm; // outer radius
	G4double w_end_sta_Al1 = 0.0 * deg; // start angle
	G4double w_end_spa_Al1 = 360 * deg; // span angle
	G4double w_end_x_Al1 = 0.0 * mm; // x location
	G4double w_end_y_Al1 = 0.0 * mm; // y location
	G4double w_end_z_Al1 = d9_z + d9_hd + 0.54 * mm + w_end_hd_Al1; // z location
	
	G4Tubs* solid_w_end_Al1 = new G4Tubs("W_end_plate_Al1", w_end_ir_Al1, w_end_or_Al1, w_end_hd_Al1, w_end_sta_Al1, w_end_spa_Al1);
	logic_w_end_Al1 = new G4LogicalVolume(solid_w_end_Al1, Alalloy, "W_end_Al1", 0, 0, 0);
	physi_w_end_Al1 = new G4PVPlacement(0, G4ThreeVector(w_end_x_Al1, w_end_y_Al1, w_end_z_Al1), logic_w_end_Al1, "W_end_Al1", logic_w, false, 0);

	// Al Insert: Shims placed behind rear tungsten shielding
	G4double w_end_d_Al2 = 1 * mm; // depth
	G4double w_end_hd_Al2 = 0.5 * w_end_d_Al2 * mm; // half depth
	G4double w_end_ir_Al2 = 0.0 * mm; // inner radius
	G4double w_end_or_Al2 = 34.75 * mm; // outer radius
	G4double w_end_sta_Al2 = 0.0 * deg; // start angle
	G4double w_end_spa_Al2 = 360 * deg; // span angle
	G4double w_end_x_Al2 = 0.0 * mm; // x location
	G4double w_end_y_Al2 = 0.0 * mm; // y location
	G4double w_end_z_Al2 = w_end_z + w_end_hd + w_end_hd_Al2; // z location
	
	G4Tubs* solid_w_end_Al2 = new G4Tubs("W_end_plate_Al2", w_end_ir_Al2, w_end_or_Al2, w_end_hd_Al2, w_end_sta_Al2, w_end_spa_Al2);
	logic_w_end_Al2 = new G4LogicalVolume(solid_w_end_Al2, Alalloy, "W_end_Al2", 0, 0, 0);
	physi_w_end_Al2 = new G4PVPlacement(0, G4ThreeVector(w_end_x_Al2, w_end_y_Al2, w_end_z_Al2), logic_w_end_Al2, "W_end_Al2", logic_w, false, 0);

	//____________________ visible attributes ____________________

    VisAtt_w = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
    VisAtt_w->SetVisibility(true);
    logic_w->SetVisAttributes(VisAtt_w);
    
	//collimator and its teeth
	VisAtt_frontcoll = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_frontcoll->SetVisibility(true);
	logic_frontcoll->SetVisAttributes(VisAtt_frontcoll);

	VisAtt_coll = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_coll->SetVisibility(true);
	logic_coll->SetVisAttributes(VisAtt_coll);

	VisAtt_coll_embed5 = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));  //embeded heavy shielding
	VisAtt_coll_embed5->SetVisibility(true);
	logic_coll_embed5->SetVisAttributes(VisAtt_coll_embed5);

	VisAtt_coll_tooth_W1 = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));  //W1
	VisAtt_coll_tooth_W1->SetVisibility(true);
	logic_coll_tooth_W1->SetVisAttributes(VisAtt_coll_tooth_W1);

	VisAtt_coll_tooth_W2 = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));  //W2
	VisAtt_coll_tooth_W2->SetVisibility(true);
	logic_coll_tooth_W2->SetVisAttributes(VisAtt_coll_tooth_W2);

	VisAtt_coll_tooth_W3 = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));  //W3
	VisAtt_coll_tooth_W3->SetVisibility(true);
	logic_coll_tooth_W3->SetVisAttributes(VisAtt_coll_tooth_W3);

	VisAtt_coll_tooth_W4 = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));  //W4
	VisAtt_coll_tooth_W4->SetVisibility(true);
	logic_coll_tooth_W4->SetVisAttributes(VisAtt_coll_tooth_W4);

	VisAtt_coll_tooth_W5 = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));  //W5
	VisAtt_coll_tooth_W5->SetVisibility(true);
	logic_coll_tooth_W5->SetVisAttributes(VisAtt_coll_tooth_W5);

	VisAtt_al_ann1 = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_al_ann1->SetVisibility(true);
	logic_al_ann1->SetVisAttributes(VisAtt_al_ann1);

	VisAtt_al_chm = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_al_chm->SetVisibility(true);
	logic_al_chm->SetVisAttributes(VisAtt_al_chm);

	VisAtt_al_end = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_al_end->SetVisibility(true);
	logic_al_end->SetVisAttributes(VisAtt_al_end);

	VisAtt_be = new G4VisAttributes(G4Colour(0.0, 0.0, 1.0));
	VisAtt_be->SetVisibility(true);
	logic_BeWin->SetVisAttributes(VisAtt_be);

	VisAtt_w_ann = new G4VisAttributes(G4Colour(0.5, 0.0, 1.0));
	VisAtt_w_ann->SetVisibility(true);
	logic_w_ann->SetVisAttributes(VisAtt_w_ann);
	logic_w_ann1->SetVisAttributes(VisAtt_w_ann);

	VisAtt_w_ann_Al1 = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_w_ann_Al1->SetVisibility(true);
	logic_w_ann_Al1->SetVisAttributes(VisAtt_w_ann_Al1);

	VisAtt_w_chm = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
	VisAtt_w_chm->SetVisibility(true);
	logic_w_chm->SetVisAttributes(VisAtt_w_chm);

	VisAtt_w_end = new G4VisAttributes(G4Colour(0.5, 0.0, 1.0));
	VisAtt_w_end->SetVisibility(true);
	logic_w_end->SetVisAttributes(VisAtt_w_end);

	VisAtt_w_end_Al1 = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_w_end_Al1->SetVisibility(true);
	logic_w_end_Al1->SetVisAttributes(VisAtt_w_end_Al1);

	VisAtt_w_end_Al2 = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
	VisAtt_w_end_Al2->SetVisibility(true);
	logic_w_end_Al2->SetVisAttributes(VisAtt_w_end_Al2);
  
	VisAtt_d1 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d1->SetVisibility(true);
	logic_d1->SetVisAttributes(VisAtt_d1);

	VisAtt_d2 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d2->SetVisibility(true);
	logic_d2->SetVisAttributes(VisAtt_d2);

	VisAtt_d3 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d3->SetVisibility(true);
	logic_d3->SetVisAttributes(VisAtt_d3);

	VisAtt_d4 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d4->SetVisibility(true);
	logic_d4->SetVisAttributes(VisAtt_d4);

	VisAtt_d5 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d5->SetVisibility(true);
	logic_d5->SetVisAttributes(VisAtt_d5);

	VisAtt_d6 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d6->SetVisibility(true);
	logic_d6->SetVisAttributes(VisAtt_d6);

	VisAtt_d7 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d7->SetVisibility(true);
	logic_d7->SetVisAttributes(VisAtt_d7);

	VisAtt_d8 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d8->SetVisibility(true);
	logic_d8->SetVisAttributes(VisAtt_d8);

	VisAtt_d9 = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
	VisAtt_d9->SetVisibility(true);
	logic_d9->SetVisAttributes(VisAtt_d9);

// designate as sensitive detectors 
sdManager = G4SDManager::GetSDMpointer();

 //BEGIN NEW CODE FOR 8 SENSITIVE DETECTORS
 //BY DLT 3 APRIL 2008
 //Revised codes for 16 sensitive detectors
 //By H. Zhao Aug 17, 2020

 //Modified CODE FOR 9 SENSITIVE DETECTORS for Whole config
 //BY Skyler Krantz 6/8/2022
 
 G4String SiliconSensDetNames[17] = {"detector_1","detector_2","detector_3","detector_4","detector_5","detector_6","detector_7","detector_8","detector_9"};

 SiSD[0] = new SetSensDet(SiliconSensDetNames[0]);
 SiSD[1] = new SetSensDet(SiliconSensDetNames[1]);
 SiSD[2] = new SetSensDet(SiliconSensDetNames[2]);
 SiSD[3] = new SetSensDet(SiliconSensDetNames[3]);
 SiSD[4] = new SetSensDet(SiliconSensDetNames[4]);
 SiSD[5] = new SetSensDet(SiliconSensDetNames[5]);
 SiSD[6] = new SetSensDet(SiliconSensDetNames[6]);
 SiSD[7] = new SetSensDet(SiliconSensDetNames[7]);
 SiSD[8] = new SetSensDet(SiliconSensDetNames[8]);

 sdManager->AddNewDetector(SiSD[0]);
 sdManager->AddNewDetector(SiSD[1]);
 sdManager->AddNewDetector(SiSD[2]);
 sdManager->AddNewDetector(SiSD[3]);
 sdManager->AddNewDetector(SiSD[4]);
 sdManager->AddNewDetector(SiSD[5]);
 sdManager->AddNewDetector(SiSD[6]);
 sdManager->AddNewDetector(SiSD[7]);
 sdManager->AddNewDetector(SiSD[8]);

 logic_d1->SetSensitiveDetector(SiSD[0]);
 logic_d2->SetSensitiveDetector(SiSD[1]);
 logic_d3->SetSensitiveDetector(SiSD[2]);     
 logic_d4->SetSensitiveDetector(SiSD[3]);
 logic_d5->SetSensitiveDetector(SiSD[4]);
 logic_d6->SetSensitiveDetector(SiSD[5]);
 logic_d7->SetSensitiveDetector(SiSD[6]);
 logic_d8->SetSensitiveDetector(SiSD[7]);
 logic_d9->SetSensitiveDetector(SiSD[8]);
	return physi_w;
}

void DetectorConstruction::SetWindowMaterial(G4String materialChoice) 
{

  if (materialChoice == "Al") {
	  logic_BeWin-> SetMaterial(windowMaterialAl);
     
     G4cout 
          << G4endl 
          << "----> The window is made of " << materialChoice << G4endl;
    G4cout 
          << G4endl 
          << "----> The window is made of " << logic_BeWin->GetMaterial() << G4endl;
    }
  if (materialChoice == "Be") {
	  logic_BeWin-> SetMaterial(windowMaterialBe);
     
     G4cout 
          << G4endl 
          << "----> The window is made of " << materialChoice << G4endl; 
  }

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    G4cout 
          << G4endl 
          << "----> The window is made of " << logic_BeWin->GetMaterial() << G4endl;
}
 

void DetectorConstruction::SetWindowDepth(G4double val) {
    if (val > 0.0) {
        windowDepth = val;
        G4RunManager::GetRunManager()->ReinitializeGeometry();
         G4cout 
          << G4endl 
          << "----> The window's depth has changed to " << windowDepth << G4endl; 
    }
    else { 
     G4cout 
          << G4endl 
          << "----> The window depth choice must be greater than 0. -Command is ignored." << G4endl;
    };

    
}
