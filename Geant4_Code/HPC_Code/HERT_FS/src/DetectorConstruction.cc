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
//TEST

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
// CADMESH //
//#define CADMESH_LEXER_VERBOSE
#include "CADMesh.hh"


DetectorConstruction::DetectorConstruction()
	:  //fDetectorMessenger(0),
	logic_w(0),logic_alHousing(0), logic_BeWin(0), logic_TaTooth5(0), logic_TaTooth4(0), logic_TaSpac4(0), logic_TaTooth3(0), logic_TaSpac3(0),
	logic_TaTooth2(0), logic_TaSpac2(0), logic_TaTooth1(0), logic_TaSpac1(0), logic_WFr1(0), logic_AlFrShim(0), logic_WFrIn1(0), logic_AlChShim1(0),
	logic_BackW(0), logic_EpoxCham(0), logic_AlBkShim1(0), logic_AlBkPlate(0), logic_AlignPin1(0), logic_AlignPin2(0), logic_AlignPin3(0),
	
	physi_alHousing(0), physi_BeWin(0), physi_TaTooth5(0), physi_TaTooth4(0), physi_TaSpac4(0), physi_TaTooth3(0), physi_TaSpac3(0),
	physi_TaTooth2(0), physi_TaSpac2(0), physi_TaTooth1(0), physi_TaSpac1(0), physi_WFr1(0), physi_AlFrShim(0), physi_AlChShim1(0),
	physi_BackW(0), physi_EpoxCham(0), physi_AlBkShim1(0), physi_AlBkPlate(0), physi_AlignPin1(0), physi_AlignPin2(0), physi_AlignPin3(0),
	physi_WFrIn1(0) 
	
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
	
	// Aluminium (6061-T6) Reference: https://www.mcmaster.com/89015k111
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



	G4Material* EpoxyBase = new G4Material("EpoxyBase", density = 1.17 * g / cm3, ncomponents = 3);
	EpoxyBase->AddElement(elH, 21); // Hydrogen
	EpoxyBase->AddElement(elO, 3); // Oxygen
	EpoxyBase->AddElement(C, 28); // Carbon

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

	// air Reference: https://apc.u-paris.fr/~franco/g4doxy4.10/html/class_materials.html
	//80   density = 1.2929e-03 * g / cm3;  // at 20 degree
	//81   G4Material * Air = new G4Material("Air", density, nel = 3,
	//	82                                    kStateGas, expTemp);
	//83   G4double ttt = 75.47 + 23.20 + 1.28;
	//84   Air->AddElement(elN, massfraction = 75.47 / ttt);
	//85   Air->AddElement(elO, massfraction = 23.20 / ttt);
	//86   Air->AddElement(elAr, massfraction = 1.28 / ttt);

	//Air?
	//G4double atomicNumber = 7.;
	//G4double massOfMole = 14*g/mole;
	//G4double vac_density = 1.5e-3*g/cm3;
	//G4double temperature = 296.*kelvin;
	//G4double pressure = 101325.*pascal;
	//G4Material* Vacuum = new G4Material("interGalactic", atomicNumber,massOfMole, vac_density, kStateGas,temperature, pressure);
		
	//Set the window material and depth
	windowMaterialBe = Be; // this is to ensure the Be and Al properties are the same as what we set and also limit the options we have for window material
	windowMaterialAl = Al;
	windowMaterial = windowMaterialBe;
	//____________________ solids ____________________

	// world volume (w)
	G4double w_fl = 100.0*cm; // full length
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
	G4double pDTheta = 180 * deg - pSTheta; // for 15 degree cap;
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

	/* W Second Layer of Shielding */
	auto mesh_WFrIn1 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/W_Fr2.obj");
	G4VSolid* solid_WFrIn1 = mesh_WFrIn1->GetSolid();
	solid_WFrIn1->SetName("solid_WFrIn1");
	G4double x_WFrIn1 = 0 * mm;
	G4double y_WFrIn1 = 0 * mm;
	G4double z_WFrIn1 = -2 * mm;

	logic_WFrIn1 = new G4LogicalVolume(solid_WFrIn1, W, "logical_WFrIn1", 0, 0, 0);
	physi_WFrIn1 = new G4PVPlacement(rotm_AlHousing, G4ThreeVector(x_WFrIn1, y_WFrIn1, z_WFrIn1), logic_WFrIn1, "physical_WFrIn1", logic_w, false, 0);

	/* Ta Last Tooth */
	auto mesh_TaTooth5 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Tooth.obj");
	G4VSolid* solid_TaTooth5 = mesh_TaTooth5->GetSolid();
	solid_TaTooth5->SetName("solid_TaTooth5");
	G4double x_TaTooth5 = 0 * mm;
	G4double y_TaTooth5 = 0 * mm;
	G4double z_TaTooth5 = d1_z-3.75*mm;
	logic_TaTooth5 = new G4LogicalVolume(solid_TaTooth5, Ta, "logical_TaTooth5", 0, 0, 0);
	physi_TaTooth5 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaTooth5, y_TaTooth5, z_TaTooth5), logic_TaTooth5, "physical_TaTooth5", logic_w, false, 0);

	/* Be Window */
	auto mesh_BeWin = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Be_Window.obj");
	G4VSolid* solid_BeWin = mesh_BeWin->GetSolid();
	solid_BeWin->SetName("solid_BeWin");
	G4double x_BeWin = 0 * mm;
	G4double y_BeWin = 0 * mm;
	G4double z_BeWin = z_TaTooth5 - 1.5 * mm;

	logic_BeWin = new G4LogicalVolume(solid_BeWin, Be, "logical_BeWin", 0, 0, 0);
	physi_BeWin = new G4PVPlacement(rotm_Y, G4ThreeVector(x_BeWin, y_BeWin, z_BeWin), logic_BeWin, "physical_BeWin", logic_w, false, 0);
	
	 /* Ta 4th Spacer */
        auto mesh_TaSpac4 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Spacer_DART_Mod.obj");
        G4VSolid* solid_TaSpac4 = mesh_TaSpac4->GetSolid();
        solid_TaSpac4->SetName("solid_TaSpac4");
        G4double x_TaSpac4 = 0 * mm;
        G4double y_TaSpac4 = 0 * mm;
        G4double z_TaSpac4 = z_BeWin-12.5 *mm;

        logic_TaSpac4 = new G4LogicalVolume(solid_TaSpac4, Ta, "logical_TaSpac4", 0, 0, 0);
        physi_TaSpac4 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaSpac4, y_TaSpac4, z_TaSpac4), logic_TaSpac4, "physical_TaSpac4", logic_w, false, 0);

	/* Ta 4th Tooth */
	auto mesh_TaTooth4 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Tooth.obj");
	G4VSolid* solid_TaTooth4 = mesh_TaTooth4->GetSolid();
	solid_TaTooth4->SetName("solid_TaTooth4");
	G4double x_TaTooth4 = 0 * mm;
	G4double y_TaTooth4 = 0 * mm;
	G4double z_TaTooth4 = z_TaSpac4 - 1.0 * mm;
	logic_TaTooth4 = new G4LogicalVolume(solid_TaTooth4, Ta, "logical_TaTooth4", 0, 0, 0);
	physi_TaTooth4 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaTooth4, y_TaTooth4, z_TaTooth4), logic_TaTooth4, "physical_TaTooth4", logic_w, false, 0);

	/* Ta 3rd Spacer */
	auto mesh_TaSpac3 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Spacer1.obj");
	G4VSolid* solid_TaSpac3 = mesh_TaSpac3->GetSolid();
	solid_TaSpac3->SetName("solid_TaSpac3");
	G4double x_TaSpac3 = 0 * mm;
	G4double y_TaSpac3 = 0 * mm;
	G4double z_TaSpac3 = z_TaTooth4 -14.00* mm;

	logic_TaSpac3 = new G4LogicalVolume(solid_TaSpac3, Ta, "logical_TaSpac3", 0, 0, 0);
	physi_TaSpac3 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaSpac3, y_TaSpac3, z_TaSpac3), logic_TaSpac3, "physical_TaSpac3", logic_w, false, 0);

	/* Ta 3rd Tooth */
	auto mesh_TaTooth3 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Tooth.obj");
	G4VSolid* solid_TaTooth3 = mesh_TaTooth3->GetSolid();
	solid_TaTooth3->SetName("solid_TaTooth3");
	G4double x_TaTooth3 = 0 * mm;
	G4double y_TaTooth3 = 0 * mm;
	G4double z_TaTooth3 = z_TaSpac3 - 1.0 * mm;
	logic_TaTooth3 = new G4LogicalVolume(solid_TaTooth3, Ta, "logical_TaTooth3", 0, 0, 0);
	physi_TaTooth3 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaTooth3, y_TaTooth3, z_TaTooth3), logic_TaTooth3, "physical_TaTooth3", logic_w, false, 0);

	/* Ta 2nd Spacer */
	auto mesh_TaSpac2 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Spacer1.obj");
	G4VSolid* solid_TaSpac2 = mesh_TaSpac2->GetSolid();
	solid_TaSpac2->SetName("solid_TaSpac2");
	G4double x_TaSpac2 = 0 * mm;
	G4double y_TaSpac2 = 0 * mm;
	G4double z_TaSpac2 = z_TaTooth3 - 14.00 * mm;

	logic_TaSpac2 = new G4LogicalVolume(solid_TaSpac2, Ta, "logical_TaSpac2", 0, 0, 0);
	physi_TaSpac2 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaSpac2, y_TaSpac2, z_TaSpac2), logic_TaSpac2, "physical_TaSpac2", logic_w, false, 0);

	/* Ta 2nd Tooth */
	auto mesh_TaTooth2 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Tooth.obj");
	G4VSolid* solid_TaTooth2 = mesh_TaTooth2->GetSolid();
	solid_TaTooth2->SetName("solid_TaTooth2");
	G4double x_TaTooth2 = 0 * mm;
	G4double y_TaTooth2 = 0 * mm;
	G4double z_TaTooth2 = z_TaSpac2 - 1.0 * mm;
	logic_TaTooth2 = new G4LogicalVolume(solid_TaTooth2, Ta, "logical_TaTooth2", 0, 0, 0);
	physi_TaTooth2 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaTooth2, y_TaTooth2, z_TaTooth2), logic_TaTooth2, "physical_TaTooth2", logic_w, false, 0);

	/* Ta 1st Spacer */
	auto mesh_TaSpac1 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Spacer1.obj");
	G4VSolid* solid_TaSpac1 = mesh_TaSpac1->GetSolid();
	solid_TaSpac1->SetName("solid_TaSpac1");
	G4double x_TaSpac1 = 0 * mm;
	G4double y_TaSpac1 = 0 * mm;
	G4double z_TaSpac1 = z_TaTooth2 - 14.00 * mm;

	logic_TaSpac1 = new G4LogicalVolume(solid_TaSpac1, Ta, "logical_TaSpac1", 0, 0, 0);
	physi_TaSpac1 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaSpac1, y_TaSpac1, z_TaSpac1), logic_TaSpac1, "physical_TaSpac1", logic_w, false, 0);

	/* Ta 1st Tooth */
	auto mesh_TaTooth1 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Ta_Tooth.obj");
	G4VSolid* solid_TaTooth1 = mesh_TaTooth1->GetSolid();
	solid_TaTooth1->SetName("solid_TaTooth1");
	G4double x_TaTooth1 = 0 * mm;
	G4double y_TaTooth1 = 0 * mm;
	G4double z_TaTooth1 = z_TaSpac1 - 1.0 * mm;
	logic_TaTooth1 = new G4LogicalVolume(solid_TaTooth1, Ta, "logical_TaTooth1", 0, 0, 0);
	physi_TaTooth1 = new G4PVPlacement(rotm_negX, G4ThreeVector(x_TaTooth1, y_TaTooth1, z_TaTooth1), logic_TaTooth1, "physical_TaTooth1", logic_w, false, 0);

	/* Al Front Shim */
	auto mesh_AlFrShim = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Al_Front_Shim.obj");
	G4VSolid* solid_AlFrShim = mesh_AlFrShim->GetSolid();
	solid_AlFrShim->SetName("solid_AlFrShim");
	G4double x_AlFrShim = 0 * mm;
	G4double y_AlFrShim = 0 * mm;
	G4double z_AlFrShim = z_WFrIn1 - 0.5 * mm;
	logic_AlFrShim = new G4LogicalVolume(solid_AlFrShim, Alalloy, "logical_AlFrShim", 0, 0, 0);
	physi_AlFrShim = new G4PVPlacement(rotm_AlHousing, G4ThreeVector(x_AlFrShim, y_AlFrShim, z_AlFrShim), logic_AlFrShim, "physical_AlFrShim", logic_w, false, 0);

	/* W Front Plate 1 */
	auto mesh_WFr1 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/W_Fr1.obj");
	G4VSolid* solid_WFr1 = mesh_WFr1->GetSolid();
	solid_WFr1->SetName("solid_WFr1");
	G4double x_WFr1 = 0 * mm;
	G4double y_WFr1 = 0 * mm;
	G4double z_WFr1 = z_AlFrShim - 3.5 * mm;
	logic_WFr1 = new G4LogicalVolume(solid_WFr1, W, "logical_WFr1", 0, 0, 0);
	physi_WFr1 = new G4PVPlacement(rotm_AlHousing, G4ThreeVector(x_WFr1, y_WFr1, z_WFr1), logic_WFr1, "physical_WFr1", logic_w, false, 0);


	/* Aluminum Housing */
	auto mesh_alHousing = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Al_Collimator_Housing.obj");
	G4VSolid* solid_alHousing = mesh_alHousing->GetSolid();
	solid_alHousing->SetName("solid_alHousing");
	G4double x_alHousing = 0 * mm;
	G4double y_alHousing = 0 * mm;
	G4double z_alHousing = -64.0 * mm;

	logic_alHousing = new G4LogicalVolume(solid_alHousing, Alalloy, "logical_alHousing", 0, 0, 0);
	physi_alHousing = new G4PVPlacement(rotm_AlHousing, G4ThreeVector(x_alHousing, y_alHousing, z_alHousing), logic_alHousing, "physical_alHousing", logic_w, false, 0);

	
	/* Epoxy Chamber*/
	auto mesh_EpoxCham = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/W_Epoxy_Chamber.obj");
	G4VSolid* solid_EpoxCham = mesh_EpoxCham->GetSolid();
	solid_EpoxCham->SetName("solid_EpoxCham");
	G4double x_EpoxCham = 0 * mm;
	G4double y_EpoxCham = 0 * mm;
	G4double z_EpoxCham = -0.5 * mm;

	G4RotationMatrix* rotm_EpoCham = new G4RotationMatrix();
	rotm_EpoCham->rotateY(90. * deg);
	rotm_EpoCham->rotateZ(90. * deg);

	logic_EpoxCham = new G4LogicalVolume(solid_EpoxCham, EpoxyTungsten, "logical_EpoxCham", 0, 0, 0);
	physi_EpoxCham = new G4PVPlacement(rotm_EpoCham, G4ThreeVector(x_EpoxCham, y_EpoxCham, z_EpoxCham), logic_EpoxCham, "physical_EpoxCham", logic_w, false, 0);


	/* Al Chamber Shim 1*/
	auto mesh_AlChShim1 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Al_Chamber_Shim_Block.obj");
	G4VSolid* solid_AlChShim1 = mesh_AlChShim1->GetSolid();
	solid_AlChShim1->SetName("solid_AlChShim1");
	G4double x_AlChShim1 = 0 * mm;
	G4double y_AlChShim1 = 0 * mm;
	G4double z_AlChShim1 = d9_z +0.54 * mm + 0.5*mm;

	G4RotationMatrix* rotm_AlCham = new G4RotationMatrix();
	rotm_AlCham->rotateX(-90. * deg);
	rotm_AlCham->rotateY(180. * deg);

	logic_AlChShim1 = new G4LogicalVolume(solid_AlChShim1, Alalloy, "logical_AlChShim1", 0, 0, 0);
	physi_AlChShim1 = new G4PVPlacement(rotm_AlCham, G4ThreeVector(x_AlChShim1, y_AlChShim1, z_AlChShim1), logic_AlChShim1, "physical_AlChShim1", logic_w, false, 0);


	/* Back W Shielding*/
	auto mesh_BackW = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/W_R1.obj");
	G4VSolid* solid_BackW = mesh_BackW->GetSolid();
	solid_BackW->SetName("solid_BackW");
	G4double x_BackW = 0 * mm;
	G4double y_BackW = 0 * mm;
	G4double z_BackW = z_EpoxCham + 5*mm +24.86 * mm;

	G4RotationMatrix* rotm_BackW = new G4RotationMatrix();
	rotm_BackW->rotateX(90. * deg);
	rotm_BackW->rotateY(90. * deg);
	logic_BackW = new G4LogicalVolume(solid_BackW, W, "logical_BackW", 0, 0, 0);
	physi_BackW = new G4PVPlacement(rotm_BackW, G4ThreeVector(x_BackW, y_BackW, z_BackW), logic_BackW, "physical_BackW", logic_w, false, 0);


	/* Al Back Shim 1*/
	auto mesh_AlBkShim1 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Al_Back_Shim_Block.obj");
	G4VSolid* solid_AlBkShim1 = mesh_AlBkShim1->GetSolid();
	solid_AlBkShim1->SetName("solid_AlBkShim1");
	G4double x_AlBkShim1 = 0 * mm;
	G4double y_AlBkShim1 = 0 * mm;
	G4double z_AlBkShim1 = z_BackW + 0 * mm;

	G4RotationMatrix* rotm_AlBackShim = new G4RotationMatrix();
	rotm_AlBackShim->rotateX(-90. * deg);
	rotm_AlBackShim->rotateY(180. * deg);
	logic_AlBkShim1 = new G4LogicalVolume(solid_AlBkShim1, Alalloy, "logical_AlBkShim1", 0, 0, 0);
	physi_AlBkShim1 = new G4PVPlacement(rotm_AlBackShim, G4ThreeVector(x_AlBkShim1, y_AlBkShim1, z_AlBkShim1), logic_AlBkShim1, "physical_AlBkShim1", logic_w, false, 0);


	/* Al Back Plate*/
	auto mesh_AlBkPlate = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Al_Back_Plate.obj");
	G4VSolid* solid_AlBkPlate = mesh_AlBkPlate->GetSolid();
	solid_AlBkPlate->SetName("solid_AlBkPlate");
	G4double x_AlBkPlate = 0 * mm;
	G4double y_AlBkPlate = 0 * mm;
	G4double z_AlBkPlate = z_AlBkShim1 + 6* mm;

	G4RotationMatrix* rotm_AlBackPlate = new G4RotationMatrix();
	rotm_AlBackPlate->rotateX(-90. * deg);
	rotm_AlBackPlate->rotateZ(180. * deg);

	logic_AlBkPlate = new G4LogicalVolume(solid_AlBkPlate, Alalloy, "logical_AlBkPlate", 0, 0, 0);
	physi_AlBkPlate = new G4PVPlacement(rotm_AlBackPlate, G4ThreeVector(x_AlBkPlate, y_AlBkPlate, z_AlBkPlate), logic_AlBkPlate, "physical_AlBkPlate", logic_w, false, 0);

	/* Stainless Steel Alignment Pin 1*/
	auto mesh_AlignPin1 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Alignment_Pin.obj");
	G4VSolid* solid_AlignPin1 = mesh_AlignPin1->GetSolid();
	solid_AlignPin1->SetName("solid_AlignPin1");
	G4double x_AlignPin1 = 25.730 * mm;
	G4double y_AlignPin1 = 0* mm;
	G4double z_AlignPin1 = -10 * mm;
	G4RotationMatrix* rotm_StAlignPin = new G4RotationMatrix();
	rotm_StAlignPin->rotateX(-90. * deg);
	//rotm_StAlignPin->rotateY(180. * deg);

	logic_AlignPin1 = new G4LogicalVolume(solid_AlignPin1, StainlessSteel, "logical_AlignPin1", 0, 0, 0);
	physi_AlignPin1 = new G4PVPlacement(rotm_StAlignPin, G4ThreeVector(x_AlignPin1, y_AlignPin1, z_AlignPin1), logic_AlignPin1, "physical_AlignPin1", logic_w, false, 0);

	/* Stainless Steel Alignment Pin 2*/
	auto mesh_AlignPin2 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Alignment_Pin.obj");
	G4VSolid* solid_AlignPin2 = mesh_AlignPin2->GetSolid();
	solid_AlignPin2->SetName("solid_AlignPin2");
	G4double x_AlignPin2 = -12.865 * mm;
	G4double y_AlignPin2 = 22.283 * mm;
	G4double z_AlignPin2 = z_AlignPin1;
	logic_AlignPin2 = new G4LogicalVolume(solid_AlignPin2, StainlessSteel, "logical_AlignPin2", 0, 0, 0);
	physi_AlignPin2 = new G4PVPlacement(rotm_StAlignPin, G4ThreeVector(x_AlignPin2, y_AlignPin2, z_AlignPin2), logic_AlignPin2, "physical_AlignPin2", logic_w, false, 0);

	/* Stainless Steel Alignment Pin 3*/
	auto mesh_AlignPin3 = CADMesh::TessellatedMesh::FromOBJ("/home/wzt0020/Geant4/HERT_Runs/src/GEANT4_HERT_Obj_Files/Alignment_Pin.obj");
	G4VSolid* solid_AlignPin3 = mesh_AlignPin3->GetSolid();
	solid_AlignPin3->SetName("solid_AlignPin3");
	G4double x_AlignPin3 = x_AlignPin2;
	G4double y_AlignPin3 = -y_AlignPin2;
	G4double z_AlignPin3 = z_AlignPin1;
	logic_AlignPin3 = new G4LogicalVolume(solid_AlignPin3, StainlessSteel, "logical_AlignPin3", 0, 0, 0);
	physi_AlignPin3 = new G4PVPlacement(rotm_StAlignPin, G4ThreeVector(x_AlignPin3, y_AlignPin3, z_AlignPin3), logic_AlignPin3, "physical_AlignPin3", logic_w, false, 0);

	//____________________ visible attributes ____________________

    VisAtt_w = new G4VisAttributes(true,G4Colour(0.0,0.0,0.0));
    logic_w->SetVisAttributes(VisAtt_w);
    
	//Detectors
	//VisAtt_detectors = new G4VisAttributes(true, G4Colour(1.0, 1.0, 0.0)); // lighter color
	VisAtt_detectors = new G4VisAttributes(true,G4Colour(0.9, 0.7, 0));
	logic_d1->SetVisAttributes(VisAtt_detectors);
	logic_d2->SetVisAttributes(VisAtt_detectors);
	logic_d3->SetVisAttributes(VisAtt_detectors);
	logic_d4->SetVisAttributes(VisAtt_detectors);
	logic_d5->SetVisAttributes(VisAtt_detectors);
	logic_d6->SetVisAttributes(VisAtt_detectors);
	logic_d7->SetVisAttributes(VisAtt_detectors);
	logic_d8->SetVisAttributes(VisAtt_detectors);
	logic_d9->SetVisAttributes(VisAtt_detectors);

    //collimator and its teeth
    VisAtt_TaColl = new G4VisAttributes(true,G4Colour(0.2, 0.8, 0.8));
    logic_TaTooth5->SetVisAttributes(VisAtt_TaColl);
	logic_TaTooth4->SetVisAttributes(VisAtt_TaColl);
	logic_TaSpac4->SetVisAttributes(VisAtt_TaColl);
	logic_TaTooth3->SetVisAttributes(VisAtt_TaColl);
	logic_TaSpac3->SetVisAttributes(VisAtt_TaColl);
	logic_TaTooth2->SetVisAttributes(VisAtt_TaColl);
	logic_TaSpac2->SetVisAttributes(VisAtt_TaColl);
	logic_TaTooth1->SetVisAttributes(VisAtt_TaColl);
	logic_TaSpac1->SetVisAttributes(VisAtt_TaColl);
  
	//Tungsten Parts
	VisAtt_tungsten = new G4VisAttributes(true,G4Colour(0.8, 0, 1));  //W1
	//Front Plates
	logic_WFr1->SetVisAttributes(VisAtt_tungsten);
	//logic_WFr2->SetVisAttributes(VisAtt_tungsten);
	//logic_WFr3->SetVisAttributes(VisAtt_tungsten);
	//logic_WFr4->SetVisAttributes(VisAtt_tungsten);
	//logic_WFr5->SetVisAttributes(VisAtt_tungsten);
	//logic_WFr6->SetVisAttributes(VisAtt_tungsten);
	//logic_WFr7->SetVisAttributes(VisAtt_tungsten);
	logic_WFrIn1->SetVisAttributes(VisAtt_tungsten);
	//logic_WFrIn2->SetVisAttributes(VisAtt_tungsten);
	//logic_WFrIn3->SetVisAttributes(VisAtt_tungsten);

	//Back Shielding
	logic_BackW->SetVisAttributes(VisAtt_tungsten);

	//Aluminum Parts
    VisAtt_al = new G4VisAttributes(true,G4Colour(0.35, 0.35, 0.35));

	//Housing
	logic_alHousing->SetVisAttributes(VisAtt_al);

	//Shims
	logic_AlFrShim->SetVisAttributes(VisAtt_al);
	logic_AlChShim1->SetVisAttributes(VisAtt_al);
	//logic_AlChShim2->SetVisAttributes(VisAtt_al);
	//logic_AlChShim3->SetVisAttributes(VisAtt_al);
	//logic_AlChShim4->SetVisAttributes(VisAtt_al);
    logic_AlBkShim1->SetVisAttributes(VisAtt_al);
	//logic_AlBkShim2->SetVisAttributes(VisAtt_al);

	//Back Plate
	logic_AlBkPlate->SetVisAttributes(VisAtt_al);
    
	//Be Window
    VisAtt_be = new G4VisAttributes(true,G4Colour(0.2, 0.3, 1));
    logic_BeWin->SetVisAttributes(VisAtt_be);
	
	//Epoxy Tungsten
    VisAtt_w_chm = new G4VisAttributes(true,G4Colour(1, 0.6, 1));
	logic_EpoxCham->SetVisAttributes(VisAtt_w_chm);

	//Stainless Steel Pins
	VisAtt_StStl = new G4VisAttributes(true,G4Colour(0.7, 0.7, 0.5));
	logic_AlignPin1->SetVisAttributes(VisAtt_StStl);
	logic_AlignPin2->SetVisAttributes(VisAtt_StStl);
	logic_AlignPin3->SetVisAttributes(VisAtt_StStl);

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
 //END NEW CODE FOR DETECTORS
	return physi_w;
}

void DetectorConstruction::SetWindowMaterial(G4String materialChoice) 
{

  if (materialChoice == "Al") {
	  logic_BeWin-> SetMaterial(windowMaterialAl);
     //windowMaterial = windowMaterialAl;
     
     G4cout 
          << G4endl 
          << "----> The window is made of " << materialChoice << G4endl;
    G4cout 
          << G4endl 
          << "----> The window is made of " << logic_BeWin->GetMaterial() << G4endl;
    }
  if (materialChoice == "Be") {
	  logic_BeWin-> SetMaterial(windowMaterialBe);
     //windowMaterial = windowMaterialBe;
     
     G4cout 
          << G4endl 
          << "----> The window is made of " << materialChoice << G4endl; 
  }

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout 
          << G4endl 
          << "----> The window is made of " << logic_BeWin->GetMaterial() << G4endl;
}
 

void DetectorConstruction::SetWindowDepth(G4double val) {
    if (val > 0.0) {
        windowDepth = val;
        //G4RunManager::GetRunManager()->GeometryHasBeenModified(); 
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
