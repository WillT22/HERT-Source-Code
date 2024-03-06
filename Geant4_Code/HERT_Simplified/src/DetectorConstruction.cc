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


// Revised by Skyler Krantz: 6/7/2023
// Modified Design for simplfied HERT-CAD version


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

	logic_w(0), logic_S1(0), logic_frontcoll(0), logic_coll(0), logic_coll_embed1(0), logic_coll_embed2(0), logic_coll_embed3(0), logic_coll_embed4(0),
	logic_coll_embed5(0), logic_coll_embed_Al1(0), logic_coll_embed_Al2(0), logic_coll_embed_Al3(0), logic_coll_embed_Al4(0),
	logic_coll_tooth_W1(0), logic_coll_tooth_W2(0), logic_coll_tooth_W3(0), logic_coll_tooth_W4(0),
	logic_coll_tooth_W5(0), logic_coll_tooth_W6(0), logic_coll_tooth_W7(0),
	logic_al_ann1(0), logic_al_chm(0), logic_al_end(0), logic_al_endann(0), logic_be(0),
	logic_w_ann(0), logic_w_ann1(0), logic_w_ann2(0), logic_w_ann3(0), logic_w_ann4(0), logic_w_ann_Al1(0), logic_w_ann_Al2(0), logic_w_ann_Al3(0),
	logic_w_chm(0), logic_w_chm1(0), logic_w_chm2(0), logic_w_chm3(0), logic_w_chm4(0), logic_w_chm5(0), logic_w_chm_Al1(0), logic_w_chm_Al2(0), logic_w_chm_Al3(0), logic_w_chm_Al4(0),
	logic_w_end(0), logic_w_end1(0), logic_w_end2(0), logic_w_end3(0), logic_w_end4(0), logic_w_end5(0), logic_w_end_Al1(0), logic_w_end_Al2(0), logic_w_end_Al3(0), logic_w_end_Al4(0),
	logic_w_endenh(0), logic_w_endenh1(0), logic_w_endenh2(0), logic_w_endenh3(0), logic_w_endenh_Al1(0), logic_w_endenh_Al2(0), logic_d1(0),  logic_d2(0), logic_d3(0), logic_d4(0),
	logic_d5(0),  logic_d6(0), logic_d7(0),  logic_d8(0), logic_d9(0), logic_fastener1(0), logic_fastener2(0), logic_fastener3(0),
    physi_w(0), physi_S1(0), physi_frontcoll(0), physi_coll(0), physi_coll_embed1(0), physi_coll_embed2(0), physi_coll_embed3(0), physi_coll_embed4(0)
    , physi_coll_embed5(0), physi_coll_embed_Al1(0), physi_coll_embed_Al2(0), physi_coll_embed_Al3(0), physi_coll_embed_Al4(0),
    physi_coll_tooth_W1(0), physi_coll_tooth_W2(0), physi_coll_tooth_W3(0), physi_coll_tooth_W4(0), 
    physi_coll_tooth_W5(0), physi_coll_tooth_W6(0), physi_coll_tooth_W7(0), 
    physi_al_ann1(0), physi_al_chm(0), physi_al_end(0), physi_al_endann(0), physi_be(0), 
    physi_w_ann(0),physi_w_ann1(0),physi_w_ann2(0),physi_w_ann3(0),physi_w_ann4(0), physi_w_ann_Al1(0),physi_w_ann_Al2(0),physi_w_ann_Al3(0),
    physi_w_chm(0),physi_w_chm1(0), physi_w_chm2(0), physi_w_chm3(0), physi_w_chm4(0), physi_w_chm5(0), physi_w_chm_Al1(0),physi_w_chm_Al2(0),physi_w_chm_Al3(0),physi_w_chm_Al4(0),    
    physi_w_end(0),physi_w_end1(0), physi_w_end2(0),physi_w_end3(0),physi_w_end4(0),physi_w_end5(0),physi_w_end_Al1(0),physi_w_end_Al2(0),physi_w_end_Al3(0),physi_w_end_Al4(0),
    physi_w_endenh(0),physi_w_endenh1(0),physi_w_endenh2(0),physi_w_endenh3(0),physi_w_endenh_Al1(0),physi_w_endenh_Al2(0), physi_d1(0), physi_d2(0), physi_d3(0), physi_d4(0), physi_d5(0), physi_d6(0), physi_d7(0), physi_d8(0),
     physi_d9(0), physi_fastener1(0),physi_fastener2(0),physi_fastener3(0)
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

	/*// air Reference: https://apc.u-paris.fr/~franco/g4doxy4.10/html/class_materials.html
	80   density = 1.2929e-03 * g / cm3;  // at 20 degree
	81   G4Material * Air = new G4Material("Air", density, nel = 3,
		82                                    kStateGas, expTemp);
	83   G4double ttt = 75.47 + 23.20 + 1.28;
	84   Air->AddElement(elN, massfraction = 75.47 / ttt);
	85   Air->AddElement(elO, massfraction = 23.20 / ttt);
	86   Air->AddElement(elAr, massfraction = 1.28 / ttt);*/

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
	G4double w_fl = 40.0*cm; // full length
	G4double w_hl = 0.5*w_fl; // half length
	
	G4Box* solid_w = new G4Box("world",w_hl,w_hl,w_hl);
	logic_w = new G4LogicalVolume(solid_w,Vacuum,"world",0,0,0);
	physi_w = new G4PVPlacement(0,G4ThreeVector(),logic_w,"world",0,false,0);
	




	/*
        // W/Al fastener1: removed for HERT1
        G4double fastener1_d=31.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
        G4double fastener1_hd=0.5*fastener1_d*mm; // half depth
        G4double fastener1_ir=0.00*mm; // inner radius
        G4double fastener1_or=1.625*mm; // outer radius | (5 mm thick Al collimator)
        G4double fastener1_sta=0.0*deg; // start angle
        G4double fastener1_spa=360*deg; // span angle
        G4double fastener1_x=0.0*mm; // x location
        G4double fastener1_y=25.73*mm; // y location
        G4double fastener1_z= 32.1*mm + fastener1_hd; // z location
        G4Tubs* solid_fastener1 = new G4Tubs("fastener1",fastener1_ir,fastener1_or,fastener1_hd,fastener1_sta,fastener1_spa);
        logic_fastener1 = new G4LogicalVolume(solid_fastener1,StainlessSteel,"fastener1",0,0,0);
        physi_fastener1 = new G4PVPlacement(0,G4ThreeVector(fastener1_x,fastener1_y,fastener1_z),logic_fastener1,"fastener1",logic_w,false,0);
    
        // W/Al fastener2
        G4double fastener2_d=31.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
        G4double fastener2_hd=0.5*fastener2_d*mm; // half depth
        G4double fastener2_ir=0.00*mm; // inner radius
        G4double fastener2_or=1.625*mm; // outer radius | (5 mm thick Al collimator)
        G4double fastener2_sta=0.0*deg; // start angle
        G4double fastener2_spa=360*deg; // span angle
        G4double fastener2_x= 22.28*mm; // x location
        G4double fastener2_y=-12.86*mm; // y location
        G4double fastener2_z= 32.1*mm + fastener2_hd; // z location
        G4Tubs* solid_fastener2 = new G4Tubs("fastener2",fastener2_ir,fastener2_or,fastener2_hd,fastener2_sta,fastener2_spa);
        logic_fastener2 = new G4LogicalVolume(solid_fastener2,StainlessSteel,"fastener2",0,0,0);
        physi_fastener2 = new G4PVPlacement(0,G4ThreeVector(fastener2_x,fastener2_y,fastener2_z),logic_fastener2,"fastener2",logic_w,false,0);
    
        // W/Al fastener3
        G4double fastener3_d=31.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
        G4double fastener3_hd=0.5*fastener3_d*mm; // half depth
        G4double fastener3_ir=0.00*mm; // inner radius
        G4double fastener3_or=1.625*mm; // outer radius | (5 mm thick Al collimator)
        G4double fastener3_sta=0.0*deg; // start angle
        G4double fastener3_spa=360*deg; // span angle
        G4double fastener3_x=-22.28*mm; // x location
        G4double fastener3_y=-12.86*mm; // y location
        G4double fastener3_z= 32.1*mm + fastener3_hd; // z location
        G4Tubs* solid_fastener3 = new G4Tubs("fastener3",fastener3_ir,fastener3_or,fastener3_hd,fastener3_sta,fastener3_spa);
        logic_fastener3 = new G4LogicalVolume(solid_fastener3,StainlessSteel,"fastener3",0,0,0);
        physi_fastener3 = new G4PVPlacement(0,G4ThreeVector(fastener3_x,fastener3_y,fastener3_z),logic_fastener3,"fastener3",logic_w,false,0);
    */

        // aluminum front collimator
        G4double frontcoll_d=1.00*mm; // depth | 40.0 - 3.5 - 5.0 mm
        G4double frontcoll_hd=0.5*frontcoll_d*mm; // half depth
        G4double frontcoll_ir=9.00*mm; // inner radius
        G4double frontcoll_or=23.5*mm; // outer radius | (5 mm thick Al collimator)
        G4double frontcoll_sta=0.0*deg; // start angle
        G4double frontcoll_spa=360*deg; // span angle
        G4double frontcoll_x=0.0*mm; // x location
        G4double frontcoll_y=0.0*mm; // y location
        G4double frontcoll_z= frontcoll_hd; // z location
        G4Tubs* solid_frontcoll = new G4Tubs("Al_frontcollimator",frontcoll_ir,frontcoll_or,frontcoll_hd,frontcoll_sta,frontcoll_spa);
        logic_frontcoll = new G4LogicalVolume(solid_frontcoll,Al,"Al_frontcollimator",0,0,0);
        physi_frontcoll = new G4PVPlacement(0,G4ThreeVector(frontcoll_x,frontcoll_y,frontcoll_z),logic_frontcoll,"Al_frontcollimator",logic_w,false,0);
    
	// aluminum collimator
	G4double coll_d=47*mm; // depth
	G4double coll_hd=0.5*coll_d*mm; // half depth
	G4double coll_ir=17.125*mm; // inner radius
	G4double coll_or=23.5*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta=0.0*deg; // start angle
	G4double coll_spa=360*deg; // span angle
	G4double coll_x=0.0*mm; // x location
	G4double coll_y=0.0*mm; // y location
	G4double coll_z= frontcoll_d+coll_hd; // z location
	G4Tubs* solid_coll = new G4Tubs("Al_collimator",coll_ir,coll_or,coll_hd,coll_sta,coll_spa);
        logic_coll = new G4LogicalVolume(solid_coll,Al,"Al_collimator",0,0,0);
        physi_coll = new G4PVPlacement(0,G4ThreeVector(coll_x,coll_y,coll_z),logic_coll,"Al_collimator",logic_w,false,0);

/*
        // heavy shielding embedded in collimator part 1
	G4double coll_d_embed1    =6.9768*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed1   =0.5*coll_d_embed1*mm; // half depth
	G4double coll_ir_embed1   =15.0*mm; // inner radius
	G4double coll_or_embed1   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed1   =0.0*deg; // start angle
	G4double coll_spa_embed1   =360*deg; // span angle
	G4double coll_x_embed1     =0.0*mm; // x location
	G4double coll_y_embed1     =0.0*mm; // y location
	G4double coll_z_embed1     =frontcoll_d+coll_hd_embed1; // z location
	G4Tubs* solid_coll_embed1  = new G4Tubs("embed_collimator1",coll_ir_embed1,coll_or_embed1,
        coll_hd_embed1,coll_sta_embed1,coll_spa_embed1);
        logic_coll_embed1        = new G4LogicalVolume(solid_coll_embed1,Ta,"embed_collimator1",0,0,0);
        physi_coll_embed1         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed1,coll_y_embed1,coll_z_embed1),
        logic_coll_embed1,"embed_collimator1",logic_w,false,0);
	// Al inserts heavy shielding embedded in collimator part 1
	G4double coll_d_embed_Al1    =0.404*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed_Al1   =0.5*coll_d_embed_Al1*mm; // half depth
	G4double coll_ir_embed_Al1   =15.0*mm; // inner radius
	G4double coll_or_embed_Al1   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed_Al1   =0.0*deg; // start angle
	G4double coll_spa_embed_Al1   =360*deg; // span angle
	G4double coll_x_embed_Al1     =0.0*mm; // x location
	G4double coll_y_embed_Al1     =0.0*mm; // y location
	G4double coll_z_embed_Al1     =frontcoll_d+coll_d_embed1+coll_hd_embed_Al1; // z location
	G4Tubs* solid_coll_embed_Al1  = new G4Tubs("embed_collimator_Al1",coll_ir_embed_Al1,coll_or_embed_Al1,
    coll_hd_embed_Al1,coll_sta_embed_Al1,coll_spa_embed_Al1);
    logic_coll_embed_Al1        = new G4LogicalVolume(solid_coll_embed_Al1,Alalloy,"embed_collimator_Al1",0,0,0);
    physi_coll_embed_Al1         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed_Al1,coll_y_embed_Al1,coll_z_embed_Al1),
    logic_coll_embed_Al1,"embed_collimator_Al1",logic_w,false,0);
        
        
 	// heavy shielding embedded in collimator part 2
	G4double coll_d_embed2    =6.9768*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed2   =0.5*coll_d_embed2*mm; // half depth
	G4double coll_ir_embed2   =15.0*mm; // inner radius
	G4double coll_or_embed2   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed2   =0.0*deg; // start angle
	G4double coll_spa_embed2   =360*deg; // span angle
	G4double coll_x_embed2     =0.0*mm; // x location
	G4double coll_y_embed2     =0.0*mm; // y location
	G4double coll_z_embed2     =frontcoll_d+coll_d_embed1+coll_d_embed_Al1+coll_hd_embed2; // z location
	G4Tubs* solid_coll_embed2  = new G4Tubs("embed_collimator2",coll_ir_embed2,coll_or_embed2,
    coll_hd_embed2,coll_sta_embed2,coll_spa_embed2);
    logic_coll_embed2        = new G4LogicalVolume(solid_coll_embed2,Ta,"embed_collimator2",0,0,0);
    physi_coll_embed2         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed2,coll_y_embed2,coll_z_embed2),
    logic_coll_embed2,"embed_collimator2",logic_w,false,0);
    // Al Insert: heavy shielding embedded in collimator part 2
    G4double coll_d_embed_Al2   =0.404*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed_Al2   =0.5*coll_d_embed_Al2*mm; // half depth
	G4double coll_ir_embed_Al2   =15.0*mm; // inner radius
	G4double coll_or_embed_Al2   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed_Al2   =0.0*deg; // start angle
	G4double coll_spa_embed_Al2   =360*deg; // span angle
	G4double coll_x_embed_Al2     =0.0*mm; // x location
	G4double coll_y_embed_Al2     =0.0*mm; // y location
	G4double coll_z_embed_Al2     =frontcoll_d+coll_d_embed1+coll_d_embed_Al1+coll_d_embed2+coll_hd_embed_Al2; // z location
	G4Tubs* solid_coll_embed_Al2  = new G4Tubs("embed_collimator_Al2",coll_ir_embed_Al2,coll_or_embed_Al2,
    coll_hd_embed_Al2,coll_sta_embed_Al2,coll_spa_embed_Al2);
    logic_coll_embed_Al2         = new G4LogicalVolume(solid_coll_embed_Al2,Alalloy,"embed_collimator_Al2",0,0,0);
    physi_coll_embed_Al2         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed_Al2,coll_y_embed_Al2,coll_z_embed_Al2),
    logic_coll_embed_Al2,"embed_collimator_Al2",logic_w,false,0);
        
    // heavy shielding embedded in collimator part 3
	G4double coll_d_embed3    =6.9768*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed3   =0.5*coll_d_embed3*mm; // half depth
	G4double coll_ir_embed3   =15.0*mm; // inner radius
	G4double coll_or_embed3   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed3   =0.0*deg; // start angle
	G4double coll_spa_embed3   =360*deg; // span angle
	G4double coll_x_embed3     =0.0*mm; // x location
	G4double coll_y_embed3     =0.0*mm; // y location
	G4double coll_z_embed3     =frontcoll_d+coll_d_embed1+coll_d_embed_Al1+coll_d_embed2+coll_d_embed_Al2+coll_hd_embed3; // z location
	G4Tubs* solid_coll_embed3  = new G4Tubs("embed_collimator3",coll_ir_embed3,coll_or_embed3,
    coll_hd_embed3,coll_sta_embed3,coll_spa_embed3);
    logic_coll_embed3        = new G4LogicalVolume(solid_coll_embed3,Ta,"embed_collimator3",0,0,0);
    physi_coll_embed3         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed3,coll_y_embed3,coll_z_embed3),
    logic_coll_embed3,"embed_collimator3",logic_w,false,0);
    // Al Insert: heavy shielding embedded in collimator part 3
    G4double coll_d_embed_Al3   =0.404*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed_Al3   =0.5*coll_d_embed_Al3*mm; // half depth
	G4double coll_ir_embed_Al3   =15.0*mm; // inner radius
	G4double coll_or_embed_Al3   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed_Al3   =0.0*deg; // start angle
	G4double coll_spa_embed_Al3   =360*deg; // span angle
	G4double coll_x_embed_Al3     =0.0*mm; // x location
	G4double coll_y_embed_Al3     =0.0*mm; // y location
	G4double coll_z_embed_Al3     =frontcoll_d+coll_d_embed1+coll_d_embed_Al1+coll_d_embed2+coll_hd_embed_Al2+coll_d_embed3+coll_hd_embed_Al3; // z location
	G4Tubs* solid_coll_embed_Al3  = new G4Tubs("embed_collimator_Al3",coll_ir_embed_Al3,coll_or_embed_Al3,
    coll_hd_embed_Al3,coll_sta_embed_Al3,coll_spa_embed_Al3);
    logic_coll_embed_Al3         = new G4LogicalVolume(solid_coll_embed_Al3,Alalloy,"embed_collimator_Al3",0,0,0);
    physi_coll_embed_Al3         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed_Al3,coll_y_embed_Al3,coll_z_embed_Al3),
    logic_coll_embed_Al3,"embed_collimator_Al3",logic_w,false,0);
              
    // heavy shielding embedded in collimator part 4
	G4double coll_d_embed4    =6.9768*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed4   =0.5*coll_d_embed4*mm; // half depth
	G4double coll_ir_embed4   =15.0*mm; // inner radius
	G4double coll_or_embed4   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed4   =0.0*deg; // start angle
	G4double coll_spa_embed4   =360*deg; // span angle
	G4double coll_x_embed4     =0.0*mm; // x location
	G4double coll_y_embed4     =0.0*mm; // y location
	G4double coll_z_embed4     =frontcoll_d+coll_d_embed1+coll_d_embed_Al1+coll_d_embed2+coll_hd_embed_Al2+coll_d_embed3+coll_d_embed_Al3+coll_hd_embed4; // z location
	G4Tubs* solid_coll_embed4  = new G4Tubs("embed_collimator4",coll_ir_embed4,coll_or_embed4,
        coll_hd_embed4,coll_sta_embed4,coll_spa_embed4);
        logic_coll_embed4        = new G4LogicalVolume(solid_coll_embed4,Ta,"embed_collimator4",0,0,0);
        physi_coll_embed4         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed4,coll_y_embed4,coll_z_embed4),
        logic_coll_embed4,"embed_collimator4",logic_w,false,0);
        // Al Insert: heavy shielding embedded in collimator part 4
        G4double coll_d_embed_Al4   =0.404*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed_Al4   =0.5*coll_d_embed_Al4*mm; // half depth
	G4double coll_ir_embed_Al4   =15.0*mm; // inner radius
	G4double coll_or_embed_Al4   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed_Al4   =0.0*deg; // start angle
	G4double coll_spa_embed_Al4   =360*deg; // span angle
	G4double coll_x_embed_Al4     =0.0*mm; // x location
	G4double coll_y_embed_Al4     =0.0*mm; // y location
	G4double coll_z_embed_Al4     =frontcoll_d+coll_d_embed1+coll_d_embed_Al1+coll_d_embed2+coll_hd_embed_Al2+coll_d_embed3+coll_d_embed_Al3+coll_d_embed4+coll_hd_embed_Al4; // z location
	G4Tubs* solid_coll_embed_Al4  = new G4Tubs("embed_collimator_Al4",coll_ir_embed_Al4,coll_or_embed_Al4,
        coll_hd_embed_Al4,coll_sta_embed_Al4,coll_spa_embed_Al4);
        logic_coll_embed_Al4        = new G4LogicalVolume(solid_coll_embed_Al4,Alalloy,"embed_collimator_Al4",0,0,0);
        physi_coll_embed_Al4        = new G4PVPlacement(0,G4ThreeVector(coll_x_embed_Al4,coll_y_embed_Al4,coll_z_embed_Al4),
        logic_coll_embed_Al4,"embed_collimator_Al4",logic_w,false,0);
                
        // heavy shielding embedded in collimator part 5
	G4double coll_d_embed5    =6.9768*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed5   =0.5*coll_d_embed5*mm; // half depth
	G4double coll_ir_embed5   =15.0*mm; // inner radius
	G4double coll_or_embed5   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed5   =0.0*deg; // start angle
	G4double coll_spa_embed5   =360*deg; // span angle
	G4double coll_x_embed5     =0.0*mm; // x location
	G4double coll_y_embed5     =0.0*mm; // y location
	G4double coll_z_embed5     =frontcoll_d+coll_d_embed1+coll_d_embed_Al1+coll_d_embed2+coll_hd_embed_Al2+coll_d_embed3+coll_d_embed_Al3+coll_d_embed4+coll_d_embed_Al4+coll_hd_embed5; // z location
	G4Tubs* solid_coll_embed5  = new G4Tubs("embed_collimator5",coll_ir_embed5,coll_or_embed5,
        coll_hd_embed5,coll_sta_embed5,coll_spa_embed5);
        logic_coll_embed5        = new G4LogicalVolume(solid_coll_embed5,Ta,"embed_collimator5",0,0,0);
        physi_coll_embed5         = new G4PVPlacement(0,G4ThreeVector(coll_x_embed5,coll_y_embed5,coll_z_embed5),
        logic_coll_embed5,"embed_collimator5",logic_w,false,0);
 */
        // heavy shielding embedded in collimator 
	G4double coll_d_embed5    =58.5*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_embed5   =0.5*coll_d_embed5*mm; // half depth
	G4double coll_ir_embed5   =15.0*mm; // inner radius
	G4double coll_or_embed5   =17.0*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_embed5   =0.0*deg; // start angle
	G4double coll_spa_embed5   =360*deg; // span angle
	G4double coll_x_embed5     =0.0*mm; // x location
	G4double coll_y_embed5     =0.0*mm; // y location
	G4double coll_z_embed5     =frontcoll_d+coll_hd_embed5; // z location
	G4Tubs* solid_coll_embed5  = new G4Tubs("embed_collimator5",coll_ir_embed5,coll_or_embed5,coll_hd_embed5,coll_sta_embed5,coll_spa_embed5);
        logic_coll_embed5      = new G4LogicalVolume(solid_coll_embed5,Ta,"embed_collimator5",0,0,0);
        physi_coll_embed5      = new G4PVPlacement(0,G4ThreeVector(coll_x_embed5,coll_y_embed5,coll_z_embed5),logic_coll_embed5,"embed_collimator5",logic_w,false,0);
 

	// collimator tooth - W1
	G4double coll_d_tooth_W1=1.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W1=0.5*coll_d_tooth_W1*mm; // half depth
	G4double coll_ir_tooth_W1=9.00*mm; // inner radius
	G4double coll_or_tooth_W1=15.00*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W1=0.0*deg; // start angle
	G4double coll_spa_tooth_W1=360*deg; // span angle
	G4double coll_x_tooth_W1=0.0*mm; // x location
	G4double coll_y_tooth_W1=0.0*mm; // y location
	G4double coll_z_tooth_W1=frontcoll_d+coll_hd_tooth_W1; // z location
	G4Tubs* solid_coll_tooth_W1 = new G4Tubs("collimator_tooth_W1",coll_ir_tooth_W1,coll_or_tooth_W1,
        coll_hd_tooth_W1,coll_sta_tooth_W1,coll_spa_tooth_W1);
        logic_coll_tooth_W1 = new G4LogicalVolume(solid_coll_tooth_W1,Ta,"collimator_tooth_W1",0,0,0);
        physi_coll_tooth_W1 = new G4PVPlacement(0,G4ThreeVector(coll_x_tooth_W1,coll_y_tooth_W1,coll_z_tooth_W1),
        logic_coll_tooth_W1,"collimator_tooth_W1",logic_w,false,0);
	// collimator tooth - W2
	G4double coll_d_tooth_W2=1.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W2=0.5*coll_d_tooth_W2*mm; // half depth
	G4double coll_ir_tooth_W2=9.00*mm; // inner radius
	G4double coll_or_tooth_W2=15.00*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W2=0.0*deg; // start angle
	G4double coll_spa_tooth_W2=360*deg; // span angle
	G4double coll_x_tooth_W2=0.0*mm; // x location
	G4double coll_y_tooth_W2=0.0*mm; // y location
	G4double coll_z_tooth_W2=frontcoll_d+15.0*mm+coll_hd_tooth_W2; // z location
	G4Tubs* solid_coll_tooth_W2 = new G4Tubs("collimator_tooth_W2",coll_ir_tooth_W2,coll_or_tooth_W2,
        coll_hd_tooth_W2,coll_sta_tooth_W2,coll_spa_tooth_W2);
        logic_coll_tooth_W2 = new G4LogicalVolume(solid_coll_tooth_W2,Ta,"collimator_tooth_W2",0,0,0);
        physi_coll_tooth_W2 = new G4PVPlacement(0,G4ThreeVector(coll_x_tooth_W2,coll_y_tooth_W2,coll_z_tooth_W2),
        logic_coll_tooth_W2,"collimator_tooth_W2",logic_w,false,0);
	// collimator tooth - W3
	G4double coll_d_tooth_W3=1.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W3=0.5*coll_d_tooth_W3*mm; // half depth
	G4double coll_ir_tooth_W3=9.00*mm; // inner radius
	G4double coll_or_tooth_W3=15.00*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W3=0.0*deg; // start angle
	G4double coll_spa_tooth_W3=360*deg; // span angle
	G4double coll_x_tooth_W3=0.0*mm; // x location
	G4double coll_y_tooth_W3=0.0*mm; // y location
	G4double coll_z_tooth_W3=frontcoll_d+30*mm+coll_hd_tooth_W3; // z location
	G4Tubs* solid_coll_tooth_W3 = new G4Tubs("collimator_tooth_W3",coll_ir_tooth_W3,coll_or_tooth_W3,
        coll_hd_tooth_W3,coll_sta_tooth_W3,coll_spa_tooth_W3);
        logic_coll_tooth_W3 = new G4LogicalVolume(solid_coll_tooth_W3,Ta,"collimator_tooth_W3",0,0,0);
        physi_coll_tooth_W3 = new G4PVPlacement(0,G4ThreeVector(coll_x_tooth_W3,coll_y_tooth_W3,coll_z_tooth_W3),
        logic_coll_tooth_W3,"collimator_tooth_W3",logic_w,false,0);
	// collimator tooth - W4
	G4double coll_d_tooth_W4=1.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W4=0.5*coll_d_tooth_W4*mm; // half depth
	G4double coll_ir_tooth_W4=9.00*mm; // inner radius
	G4double coll_or_tooth_W4=15.00*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W4=0.0*deg; // start angle
	G4double coll_spa_tooth_W4=360*deg; // span angle
	G4double coll_x_tooth_W4=0.0*mm; // x location
	G4double coll_y_tooth_W4=0.0*mm; // y location
	G4double coll_z_tooth_W4=frontcoll_d+45*mm+coll_hd_tooth_W4; // z location
	G4Tubs* solid_coll_tooth_W4 = new G4Tubs("collimator_tooth_W4",coll_ir_tooth_W4,coll_or_tooth_W4,
                            coll_hd_tooth_W4,coll_sta_tooth_W4,coll_spa_tooth_W4);
        logic_coll_tooth_W4 = new G4LogicalVolume(solid_coll_tooth_W4,Ta,"collimator_tooth_W4",0,0,0);
        physi_coll_tooth_W4 = new G4PVPlacement(0,G4ThreeVector(coll_x_tooth_W4,coll_y_tooth_W4,coll_z_tooth_W4),
        logic_coll_tooth_W4,"collimator_tooth_W4",logic_w,false,0);
	// collimator tooth - W5
	G4double coll_d_tooth_W5=1.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W5=0.5*coll_d_tooth_W5*mm; // half depth
	G4double coll_ir_tooth_W5=9.00*mm; // inner radius
	G4double coll_or_tooth_W5=17.00*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W5=0.0*deg; // start angle
	G4double coll_spa_tooth_W5=360*deg; // span angle
	G4double coll_x_tooth_W5=0.0*mm; // x location
	G4double coll_y_tooth_W5=0.0*mm; // y location
	G4double coll_z_tooth_W5=frontcoll_d+60*mm+coll_hd_tooth_W5; // z location
	G4Tubs* solid_coll_tooth_W5 = new G4Tubs("collimator_tooth_W5",coll_ir_tooth_W5,coll_or_tooth_W5,
        coll_hd_tooth_W5,coll_sta_tooth_W5,coll_spa_tooth_W5);
        logic_coll_tooth_W5 = new G4LogicalVolume(solid_coll_tooth_W5,Ta,"collimator_tooth_W5",0,0,0);
        physi_coll_tooth_W5 = new G4PVPlacement(0,G4ThreeVector(coll_x_tooth_W5,coll_y_tooth_W5,coll_z_tooth_W5),
        logic_coll_tooth_W5,"collimator_tooth_W5",logic_w,false,0);


	/* collimator tooth - W6
	G4double coll_d_tooth_W6=1.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W6=0.5*coll_d_tooth_W6*mm; // half depth
	G4double coll_ir_tooth_W6=10.00*mm; // inner radius
	G4double coll_or_tooth_W6=15.00*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W6=0.0*deg; // start angle
	G4double coll_spa_tooth_W6=360*deg; // span angle
	G4double coll_x_tooth_W6=0.0*mm; // x location
	G4double coll_y_tooth_W6=0.0*mm; // y location
	G4double coll_z_tooth_W6=frontcoll_d+67.5*mm+coll_hd_tooth_W6; // z location
	G4Tubs* solid_coll_tooth_W6 = new G4Tubs("collimator_tooth_W6",coll_ir_tooth_W6,coll_or_tooth_W6,
        coll_hd_tooth_W6,coll_sta_tooth_W6,coll_spa_tooth_W6);
        logic_coll_tooth_W6 = new G4LogicalVolume(solid_coll_tooth_W6,Ta,"collimator_tooth_W6",0,0,0);
        physi_coll_tooth_W6 = new G4PVPlacement(0,G4ThreeVector(coll_x_tooth_W6,coll_y_tooth_W6,coll_z_tooth_W6),
        logic_coll_tooth_W6,"collimator_tooth_W6",logic_w,false,0);
	// collimator tooth - W7
	G4double coll_d_tooth_W7=1.0*mm; // depth | 40.0 - 3.5 - 5.0 mm
	G4double coll_hd_tooth_W7=0.5*coll_d_tooth_W7*mm; // half depth
	G4double coll_ir_tooth_W7=10.00*mm; // inner radius
	G4double coll_or_tooth_W7=15.00*mm; // outer radius | (5 mm thick Al collimator)
	G4double coll_sta_tooth_W7=0.0*deg; // start angle
	G4double coll_spa_tooth_W7=360*deg; // span angle
	G4double coll_x_tooth_W7=0.0*mm; // x location
	G4double coll_y_tooth_W7=0.0*mm; // y location
	G4double coll_z_tooth_W7=frontcoll_d+5.5*mm+coll_hd_tooth_W7; // z location
	G4Tubs* solid_coll_tooth_W7 = new G4Tubs("collimator_tooth_W7",coll_ir_tooth_W7,coll_or_tooth_W7,
        coll_hd_tooth_W7,coll_sta_tooth_W7,coll_spa_tooth_W7);
        logic_coll_tooth_W7 = new G4LogicalVolume(solid_coll_tooth_W7,Ta,"collimator_tooth_W7",0,0,0);
        physi_coll_tooth_W7 = new G4PVPlacement(0,G4ThreeVector(coll_x_tooth_W7,coll_y_tooth_W7,coll_z_tooth_W7),
        logic_coll_tooth_W7,"collimator_tooth_W7",logic_w,false,0);        */

	// aluminum front annulus 1 (larger of the two)
	G4double al_ann1_d=10*mm; // depth
	G4double al_ann1_hd=0.5*al_ann1_d*mm; // half depth
	G4double al_ann1_ir=17.0*mm; // inner radius
	G4double al_ann1_or=40*mm; // outer radius
	G4double al_ann1_sta=0.0*deg; // start angle
	G4double al_ann1_spa=360*deg; // span angle
	G4double al_ann1_x=0.0*mm; // x location
	G4double al_ann1_y=0.0*mm; // y location
	G4double al_ann1_z=frontcoll_d+coll_d+al_ann1_hd; // z location
	G4Tubs* solid_al_ann1 = new G4Tubs("Al_annulus_1",al_ann1_ir,al_ann1_or,al_ann1_hd,al_ann1_sta,al_ann1_spa);
    //    G4VSolid* solid_al_ann1_subfastener1 = new G4SubtractionSolid("Al_annulus_subfast1",solid_al_ann1,solid_fastener1);
    //    G4VSolid* solid_al_ann1_subfastener2 = new G4SubtractionSolid("Al_annulus_subfast2",solid_al_ann1_subfastener1,solid_fastener2);
    //    G4VSolid* solid_al_ann1_subfastener3 = new G4SubtractionSolid("Al_annulus_subfast3",solid_al_ann1_subfastener2,solid_fastener3);

        logic_al_ann1 = new G4LogicalVolume(solid_al_ann1,Al,"Al_annulus_1",0,0,0);
        physi_al_ann1 = new G4PVPlacement(0,G4ThreeVector(al_ann1_x,al_ann1_y,al_ann1_z),logic_al_ann1,"Al_annulus_1",logic_w,false,0);

	// aluminum front annulus 2: holds beryllium disc (smaller of the two)
	//G4double al_ann2_d=0.5*mm; // depth
	//G4double al_ann2_hd=0.5*al_ann2_d*mm; // half depth
//	G4double al_ann2_ir=16.0*mm; // inner radius
//	G4double al_ann2_or=30.5625*mm; // outer radius
//	G4double al_ann2_sta=0.0*deg; // start angle
//	G4double al_ann2_spa=360*deg; // span angle
//	G4double al_ann2_x=0.0*mm; // x location
//	G4double al_ann2_y=0.0*mm; // y location
//	G4double al_ann2_z=coll_d+al_ann1_d+al_ann2_hd; // z location
//	G4Tubs* solid_al_ann2 = new G4Tubs("Al_annulus_2",al_ann2_ir,al_ann2_or,al_ann2_hd,al_ann2_sta,al_ann2_spa);
 //   logic_al_ann2 = new G4LogicalVolume(solid_al_ann2,Al,"Al_annulus_2",0,0,0);
  //  physi_al_ann2 = new G4PVPlacement(0,G4ThreeVector(al_ann2_x,al_ann2_y,al_ann2_z),logic_al_ann2,"Al_annulus_2",logic_w,false,0);

	// aluminum chamber
	G4double al_chm_d=36.5*mm; // depth: revised for HERT1
	G4double al_chm_hd=0.5*al_chm_d*mm; // half depth
	G4double al_chm_ir=35*mm; // inner radius | 30.5625 - 5.0 - 3.5 = 22.0625
	G4double al_chm_or=40*mm; // outer radius
	G4double al_chm_sta=0.0*deg; // start angle
	G4double al_chm_spa=360*deg; // span angle
	G4double al_chm_x=0.0*mm; // x location
	G4double al_chm_y=0.0*mm; // y location
	G4double al_chm_z=frontcoll_d+coll_d+al_ann1_d+al_chm_hd; // z location. Add 0.75 for width of delrin lip
	G4Tubs* solid_al_chm = new G4Tubs("Al_chamber",al_chm_ir,al_chm_or,al_chm_hd,al_chm_sta,al_chm_spa);
    logic_al_chm = new G4LogicalVolume(solid_al_chm,Al,"Al_chamber",0,0,0);
    physi_al_chm = new G4PVPlacement(0,G4ThreeVector(al_chm_x,al_chm_y,al_chm_z),logic_al_chm,"Al_chamber",logic_w,false,0);

	// aluminum end annulus: removed for HERT1
	//G4double al_endann_d=2.0*mm; // depth
	//G4double al_endann_hd=0.5*al_endann_d*mm; // half depth
	//G4double al_endann_ir=30.2*mm; // inner radius
	//G4double al_endann_or=38.7*mm; // outer radius
	//G4double al_endann_sta=0.0*deg; // start angle
	//G4double al_endann_spa=360*deg; // span angle
	//G4double al_endann_x=0.0*mm; // x location
	//G4double al_endann_y=0.0*mm; // y location
	//G4double al_endann_z=frontcoll_d+coll_d+al_ann1_d+al_chm_d+al_endann_hd; // z location
	//G4Tubs* solid_al_endann = new G4Tubs("Al_endann_plate",al_endann_ir,al_endann_or,
        //al_endann_hd,al_endann_sta,al_endann_spa);
        //logic_al_endann = new G4LogicalVolume(solid_al_endann,Al,"Al_endann_plate",0,0,0);
        //physi_al_endann = new G4PVPlacement(0,G4ThreeVector(al_endann_x,al_endann_y,al_endann_z),
        //logic_al_endann,"Al_endann_plate",logic_w,false,0);
	
	// aluminum end plate
	G4double al_end_d=5.0*mm; // depth
	G4double al_end_hd=0.5*al_end_d*mm; // half depth
	G4double al_end_ir=0.0*mm; // inner radius
	G4double al_end_or=40*mm; // outer radius
	G4double al_end_sta=0.0*deg; // start angle
	G4double al_end_spa=360*deg; // span angle
	G4double al_end_x=0.0*mm; // x location
	G4double al_end_y=0.0*mm; // y location
	G4double al_end_z=frontcoll_d+coll_d+al_ann1_d+al_chm_d+al_end_hd; // z location
	G4Tubs* solid_al_end = new G4Tubs("Al_end_plate",al_end_ir,al_end_or,al_end_hd,al_end_sta,al_end_spa);
    //G4VSolid* solid_al_end_subfastener1 = new G4SubtractionSolid("Al_end_plate_subfast1",solid_al_end,solid_fastener1);
    //G4VSolid* solid_al_end_subfastener2 = new G4SubtractionSolid("Al_end_plate_subfast2",solid_al_end_subfastener1,solid_fastener2);
	//G4VSolid* solid_al_end_subfastener3 = new G4SubtractionSolid("Al_end_plate_subfast3", solid_al_end_subfastener2, solid_fastener3);
    logic_al_end = new G4LogicalVolume(solid_al_end,Al,"Al_end_plate",0,0,0);
    physi_al_end = new G4PVPlacement(0,G4ThreeVector(al_end_x,al_end_y,al_end_z),logic_al_end,"Al_end_plate",logic_w,false,0);

	// beryllium disc
	//G4double be_d=0.15*mm; // depth
	G4double be_d=windowDepth; // depth
	G4double be_hd=0.5*windowDepth*mm; // half depth
	G4double be_ir=0.0*mm; // inner radius
	G4double be_or=17.0*mm; // outer radius
	G4double be_sta=0.0*deg; // start angle
	G4double be_spa=360*deg; // span angle
	G4double be_x=0.0*mm; // x location
	G4double be_y=0.0*mm; // y location
	G4double be_z=frontcoll_d+58.5*mm+0.5*windowDepth; // z location
	
	G4Tubs* solid_be = new G4Tubs("FOV_disc",be_ir,be_or,0.5*windowDepth,be_sta,be_spa);
    logic_be = new G4LogicalVolume(solid_be,windowMaterial,"FOV_disc",0,0,0);
    physi_be = new G4PVPlacement(0,G4ThreeVector(be_x,be_y,be_z),logic_be,"FOV_disc",logic_w,false,0);

    G4cout << "----> The window's z is " << be_z << G4endl; 
    G4cout << "----> The window's material is " << logic_be->GetMaterial() << G4endl;  
	G4cout << "----> The window's depth is " << windowDepth << "mm" << G4endl;

	//// tungsten front annulus Part 1
	G4double w_ann_d=3.5*mm; // depth
	G4double w_ann_hd=0.5*w_ann_d*mm; // half depth
	G4double w_ann_ir=17.125*mm; // inner radius
	G4double w_ann_or=34.75*mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta=0.0*deg; // start angle
	G4double w_ann_spa=360*deg; // span angle
	G4double w_ann_x=0.0*mm; // x location
	G4double w_ann_y=0.0*mm; // y location
	G4double w_ann_z=frontcoll_d+coll_d+al_ann1_d+w_ann_hd; // z location
	G4Tubs* solid_w_ann = new G4Tubs("W_annulus",w_ann_ir,w_ann_or,w_ann_hd,w_ann_sta,w_ann_spa);
 //       G4VSolid* solid_w_ann_subfastener1 = new G4SubtractionSolid("W_annulus_subfast1",solid_w_ann,solid_fastener1);
 //       G4VSolid* solid_w_ann_subfastener2 = new G4SubtractionSolid("W_annulus_subfast2",solid_w_ann_subfastener1,solid_fastener2);
 //       G4VSolid* solid_w_ann_subfastener3 = new G4SubtractionSolid("W_annulus_subfast3",solid_w_ann_subfastener2,solid_fastener3);
        logic_w_ann = new G4LogicalVolume(solid_w_ann, W,"W_annulus",0,0,0);
        physi_w_ann = new G4PVPlacement(0,G4ThreeVector(w_ann_x,w_ann_y,w_ann_z),logic_w_ann,"W_annulus",logic_w,false,0);
    
   
       	// Al Insert: Between tungsten part 1 and Part 2 front annulus Part 1
	G4double w_ann_d_Al1=0.5*mm; // depth
	G4double w_ann_hd_Al1=0.5*w_ann_d_Al1*mm; // half depth
	G4double w_ann_ir_Al1=17.125*mm; // inner radius
	G4double w_ann_or_Al1=34.75*mm; // outer radius | 
	G4double w_ann_sta_Al1=0.0*deg; // start angle
	G4double w_ann_spa_Al1=360*deg; // span angle
	G4double w_ann_x_Al1=0.0*mm; // x location
	G4double w_ann_y_Al1=0.0*mm; // y location
	G4double w_ann_z_Al1=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_ann_hd_Al1; // z location
	G4Tubs* solid_w_ann_Al1 = new G4Tubs("W_annulus_Al1",w_ann_ir_Al1,w_ann_or_Al1,w_ann_hd_Al1,w_ann_sta_Al1,w_ann_spa_Al1);
        //G4VSolid* solid_w_ann_Al1_subfastener1 = new G4SubtractionSolid("W_annulus_Al1_subfast1",solid_w_ann_Al1,solid_fastener1);
        //G4VSolid* solid_w_ann_Al1_subfastener2 = new G4SubtractionSolid("W_annulus_Al1_subfast2",solid_w_ann_Al1_subfastener1,solid_fastener2);
        //G4VSolid* solid_w_ann_Al1_subfastener3 = new G4SubtractionSolid("W_annulus_Al1_subfast3",solid_w_ann_Al1_subfastener2,solid_fastener3);
        logic_w_ann_Al1 = new G4LogicalVolume(solid_w_ann_Al1,Alalloy,"W_annulus_Al1",0,0,0);
        physi_w_ann_Al1 = new G4PVPlacement(0,G4ThreeVector(w_ann_x_Al1,w_ann_y_Al1,w_ann_z_Al1),logic_w_ann_Al1,"W_annulus_Al1",logic_w,false,0);
        
		// tungsten front annulus Part 2
		G4double w_ann_d1 = 1.5 * mm; // depth
		G4double w_ann_hd1 = 0.5 * w_ann_d1 * mm; // half depth
		G4double w_ann_ir1 = 11.5 * mm; // inner radius
		G4double w_ann_or1 = 34.75 * mm; // outer radius | 30.5625 - 5.0
		G4double w_ann_sta1 = 0.0 * deg; // start angle
		G4double w_ann_spa1 = 360 * deg; // span angle
		G4double w_ann_x1 = 0.0 * mm; // x location
		G4double w_ann_y1 = 0.0 * mm; // y location
		G4double w_ann_z1 = frontcoll_d + coll_d + al_ann1_d + w_ann_d + w_ann_d_Al1+w_ann_hd1; // z location
		G4Tubs* solid_w_ann1 = new G4Tubs("W_annulus1", w_ann_ir1, w_ann_or1, w_ann_hd1, w_ann_sta1, w_ann_spa1);
		//G4VSolid* solid_w_ann1_subfastener1 = new G4SubtractionSolid("W_annulus1_subfast1",solid_w_ann1,solid_fastener1);
		//G4VSolid* solid_w_ann1_subfastener2 = new G4SubtractionSolid("W_annulus1_subfast2",solid_w_ann1_subfastener1,solid_fastener2);
		//G4VSolid* solid_w_ann1_subfastener3 = new G4SubtractionSolid("W_annulus1_subfast3",solid_w_ann1_subfastener2,solid_fastener3);
		logic_w_ann1 = new G4LogicalVolume(solid_w_ann1, W, "W_annulus1", 0, 0, 0);
		physi_w_ann1 = new G4PVPlacement(0, G4ThreeVector(w_ann_x1, w_ann_y1, w_ann_z1), logic_w_ann1, "W_annulus1", logic_w, false, 0);

		/*
        // tungsten front annulus Part 2
	G4double w_ann_d2=0.572*mm; // depth
	G4double w_ann_hd2=0.5*w_ann_d2*mm; // half depth
	G4double w_ann_ir2=15.0*mm; // inner radius
	G4double w_ann_or2=33.7*mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta2=0.0*deg; // start angle
	G4double w_ann_spa2=360*deg; // span angle
	G4double w_ann_x2=0.0*mm; // x location
	G4double w_ann_y2=0.0*mm; // y location
	G4double w_ann_z2=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_hd2; // z location
	G4Tubs* solid_w_ann2 = new G4Tubs("W_annulus2",w_ann_ir2,w_ann_or2,w_ann_hd2,w_ann_sta2,w_ann_spa2);
        G4VSolid* solid_w_ann2_subfastener1 = new G4SubtractionSolid("W_annulus2_subfast1",solid_w_ann2,solid_fastener1);
        G4VSolid* solid_w_ann2_subfastener2 = new G4SubtractionSolid("W_annulus2_subfast2",solid_w_ann2_subfastener1,solid_fastener2);
        G4VSolid* solid_w_ann2_subfastener3 = new G4SubtractionSolid("W_annulus2_subfast3",solid_w_ann2_subfastener2,solid_fastener3);
        logic_w_ann2 = new G4LogicalVolume(solid_w_ann2_subfastener3,W,"W_annulus2_subfast3",0,0,0);
        physi_w_ann2 = new G4PVPlacement(0,G4ThreeVector(w_ann_x2,w_ann_y2,w_ann_z2),logic_w_ann2,"W_annulus2_subfast3",logic_w,false,0);
        
       	// Al Insert: tungsten front annulus Part 2
	G4double w_ann_d_Al2=0.404*mm; // depth
	G4double w_ann_hd_Al2=0.5*w_ann_d_Al2*mm; // half depth
	G4double w_ann_ir_Al2=15.0*mm; // inner radius
	G4double w_ann_or_Al2=33.7*mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta_Al2=0.0*deg; // start angle
	G4double w_ann_spa_Al2=360*deg; // span angle
	G4double w_ann_x_Al2=0.0*mm; // x location
	G4double w_ann_y_Al2=0.0*mm; // y location
	G4double w_ann_z_Al2=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_hd_Al2; // z location
	G4Tubs* solid_w_ann_Al2 = new G4Tubs("W_annulus_Al2",w_ann_ir_Al2,w_ann_or_Al2,w_ann_hd_Al2,w_ann_sta_Al2,w_ann_spa_Al2);
        G4VSolid* solid_w_ann_Al2_subfastener1 = new G4SubtractionSolid("W_annulus_Al2_subfast1",solid_w_ann_Al2,solid_fastener1);
        G4VSolid* solid_w_ann_Al2_subfastener2 = new G4SubtractionSolid("W_annulus_Al2_subfast2",solid_w_ann_Al2_subfastener1,solid_fastener2);
        G4VSolid* solid_w_ann_Al2_subfastener3 = new G4SubtractionSolid("W_annulus_Al2_subfast3",solid_w_ann_Al2_subfastener2,solid_fastener3);
        logic_w_ann_Al2 = new G4LogicalVolume(solid_w_ann_Al2_subfastener3,Alalloy,"W_annulus_Al2_subfast3",0,0,0);
        physi_w_ann_Al2 = new G4PVPlacement(0,G4ThreeVector(w_ann_x_Al2,w_ann_y_Al2,w_ann_z_Al2),logic_w_ann_Al2,"W_annulus_Al2_subfast3",logic_w,false,0);

        // tungsten front annulus Part 3
	G4double w_ann_d3=0.572*mm; // depth
	G4double w_ann_hd3=0.5*w_ann_d3*mm; // half depth
	G4double w_ann_ir3=15.0*mm; // inner radius
	G4double w_ann_or3=33.7*mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta3=0.0*deg; // start angle
	G4double w_ann_spa3=360*deg; // span angle
	G4double w_ann_x3=0.0*mm; // x location
	G4double w_ann_y3=0.0*mm; // y location
	G4double w_ann_z3=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_hd3; // z location
	G4Tubs* solid_w_ann3 = new G4Tubs("W_annulus3",w_ann_ir3,w_ann_or3,w_ann_hd3,w_ann_sta3,w_ann_spa3);
        G4VSolid* solid_w_ann3_subfastener1 = new G4SubtractionSolid("W_annulus3_subfast1",solid_w_ann3,solid_fastener1);
        G4VSolid* solid_w_ann3_subfastener2 = new G4SubtractionSolid("W_annulus3_subfast2",solid_w_ann3_subfastener1,solid_fastener2);
        G4VSolid* solid_w_ann3_subfastener3 = new G4SubtractionSolid("W_annulus3_subfast3",solid_w_ann3_subfastener2,solid_fastener3);
        logic_w_ann3 = new G4LogicalVolume(solid_w_ann3_subfastener3,W,"W_annulus3_subfast3",0,0,0);
        physi_w_ann3 = new G4PVPlacement(0,G4ThreeVector(w_ann_x3,w_ann_y3,w_ann_z3),logic_w_ann3,"W_annulus3_subfast3",logic_w,false,0);
        
       	// Al Insert: tungsten front annulus Part 3
	G4double w_ann_d_Al3=0.404*mm; // depth
	G4double w_ann_hd_Al3=0.5*w_ann_d_Al3*mm; // half depth
	G4double w_ann_ir_Al3=15.0*mm; // inner radius
	G4double w_ann_or_Al3=33.7*mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta_Al3=0.0*deg; // start angle
	G4double w_ann_spa_Al3=360*deg; // span angle
	G4double w_ann_x_Al3=0.0*mm; // x location
	G4double w_ann_y_Al3=0.0*mm; // y location
	G4double w_ann_z_Al3=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_hd_Al3; // z location
	G4Tubs* solid_w_ann_Al3 = new G4Tubs("W_annulus_Al3",w_ann_ir_Al3,w_ann_or_Al3,w_ann_hd_Al3,w_ann_sta_Al3,w_ann_spa_Al3);
        G4VSolid* solid_w_ann_Al3_subfastener1 = new G4SubtractionSolid("W_annulus_Al3_subfast1",solid_w_ann_Al3,solid_fastener1);
        G4VSolid* solid_w_ann_Al3_subfastener2 = new G4SubtractionSolid("W_annulus_Al3_subfast2",solid_w_ann_Al3_subfastener1,solid_fastener2);
        G4VSolid* solid_w_ann_Al3_subfastener3 = new G4SubtractionSolid("W_annulus_Al3_subfast3",solid_w_ann_Al3_subfastener2,solid_fastener3);
        logic_w_ann_Al3 = new G4LogicalVolume(solid_w_ann_Al3_subfastener3,Alalloy,"W_annulus_Al3_subfast3",0,0,0);
        physi_w_ann_Al3 = new G4PVPlacement(0,G4ThreeVector(w_ann_x_Al3,w_ann_y_Al3,w_ann_z_Al3),logic_w_ann_Al3,"W_annulus_Al3_subfast3",logic_w,false,0);

        // tungsten front annulus Part 4
	G4double w_ann_d4=0.572*mm; // depth
	G4double w_ann_hd4=0.5*w_ann_d4*mm; // half depth
	G4double w_ann_ir4=15.0*mm; // inner radius
	G4double w_ann_or4=33.7*mm; // outer radius | 30.5625 - 5.0
	G4double w_ann_sta4=0.0*deg; // start angle
	G4double w_ann_spa4=360*deg; // span angle
	G4double w_ann_x4=0.0*mm; // x location
	G4double w_ann_y4=0.0*mm; // y location
	G4double w_ann_z4=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_hd4; // z location
	G4Tubs* solid_w_ann4 = new G4Tubs("W_annulus4",w_ann_ir4,w_ann_or4,w_ann_hd4,w_ann_sta4,w_ann_spa4);
        G4VSolid* solid_w_ann4_subfastener1 = new G4SubtractionSolid("W_annulus4_subfast1",solid_w_ann4,solid_fastener1);
        G4VSolid* solid_w_ann4_subfastener2 = new G4SubtractionSolid("W_annulus4_subfast2",solid_w_ann4_subfastener1,solid_fastener2);
        G4VSolid* solid_w_ann4_subfastener3 = new G4SubtractionSolid("W_annulus4_subfast3",solid_w_ann4_subfastener2,solid_fastener3);
        logic_w_ann4 = new G4LogicalVolume(solid_w_ann4_subfastener3,W,"W_annulus4_subfast3",0,0,0);
        physi_w_ann4 = new G4PVPlacement(0,G4ThreeVector(w_ann_x4,w_ann_y4,w_ann_z4),logic_w_ann4,"W_annulus4_subfast3",logic_w,false,0);
    
	*/


	// tungsten chamber
	G4double w_chm_d=24.86*mm; // depth
	G4double w_chm_hd=0.5*w_chm_d*mm; // half depth
	G4double w_chm_ir=30.75*mm; // inner radius
	G4double w_chm_or=34.75*mm; // outer radius
	G4double w_chm_sta=0.0*deg; // start angle
	G4double w_chm_spa=360*deg; // span angle
	G4double w_chm_x=0.0*mm; // x location
	G4double w_chm_y=0.0*mm; // y location
	G4double w_chm_z= frontcoll_d + coll_d + al_ann1_d + w_ann_d + w_ann_d_Al1 + w_ann_d1+w_chm_hd; // z location
	G4Tubs* solid_w_chm = new G4Tubs("W_chamber",w_chm_ir,w_chm_or,w_chm_hd,w_chm_sta,w_chm_spa);
        logic_w_chm = new G4LogicalVolume(solid_w_chm,EpoxyTungsten,"W_chamber",0,0,0);
        physi_w_chm = new G4PVPlacement(0,G4ThreeVector(w_chm_x,w_chm_y,w_chm_z),logic_w_chm,"W_chamber",logic_w,false,0);
 
/*
        //  tungsten chamber - Part 1
	G4double w_chm_d1=1.8768*mm; // depth
	G4double w_chm_hd1=0.5*w_chm_d1*mm; // half depth
	G4double w_chm_ir1=30.2*mm; // inner radius
	G4double w_chm_or1=33.7*mm; // outer radius
	G4double w_chm_sta1=0.0*deg; // start angle
	G4double w_chm_spa1=360*deg; // span angle
	G4double w_chm_x1=0.0*mm; // x location
	G4double w_chm_y1=0.0*mm; // y location
	G4double w_chm_z1=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_hd1; // z location
	G4Tubs* solid_w_chm1 = new G4Tubs("W_chamber1",w_chm_ir1,w_chm_or1,w_chm_hd1,w_chm_sta1,w_chm_spa1);
        logic_w_chm1 = new G4LogicalVolume(solid_w_chm1,W,"W_chamber1",0,0,0);
        physi_w_chm1 = new G4PVPlacement(0,G4ThreeVector(w_chm_x1,w_chm_y1,w_chm_z1),logic_w_chm1,"W_chamber1",logic_w,false,0);
        
        // Al insert#1: tungsten chamber
	G4double w_chm_d_Al1=0.404*mm; // depth
	G4double w_chm_hd_Al1=0.5*w_chm_d_Al1*mm; // half depth
	G4double w_chm_ir_Al1=30.2*mm; // inner radius
	G4double w_chm_or_Al1=33.7*mm; // outer radius
	G4double w_chm_sta_Al1=0.0*deg; // start angle
	G4double w_chm_spa_Al1=360*deg; // span angle
	G4double w_chm_x_Al1=0.0*mm; // x location
	G4double w_chm_y_Al1=0.0*mm; // y location
	G4double w_chm_z_Al1=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_hd_Al1; // z location
	G4Tubs* solid_w_chm_Al1 = new G4Tubs("W_chamber_Al1",w_chm_ir_Al1,w_chm_or_Al1,w_chm_hd_Al1,w_chm_sta_Al1,w_chm_spa_Al1);
        logic_w_chm_Al1 = new G4LogicalVolume(solid_w_chm_Al1,W,"W_chamber_Al1",0,0,0);
        physi_w_chm_Al1 = new G4PVPlacement(0,G4ThreeVector(w_chm_x_Al1,w_chm_y_Al1,w_chm_z_Al1),logic_w_chm_Al1,"W_chamber_Al1",logic_w,false,0);
        
        // tungsten chamber Part 2
	G4double w_chm_d2=1.8768*mm; // depth
	G4double w_chm_hd2=0.5*w_chm_d2*mm; // half depth
	G4double w_chm_ir2=30.2*mm; // inner radius
	G4double w_chm_or2=33.7*mm; // outer radius
	G4double w_chm_sta2=0.0*deg; // start angle
	G4double w_chm_spa2=360*deg; // span angle
	G4double w_chm_x2=0.0*mm; // x location
	G4double w_chm_y2=0.0*mm; // y location
	G4double w_chm_z2=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_d_Al1+w_chm_hd2; // z location
	G4Tubs* solid_w_chm2 = new G4Tubs("W_chamber2",w_chm_ir2,w_chm_or2,w_chm_hd2,w_chm_sta2,w_chm_spa2);
        logic_w_chm2 = new G4LogicalVolume(solid_w_chm2,W,"W_chamber2",0,0,0);
        physi_w_chm2 = new G4PVPlacement(0,G4ThreeVector(w_chm_x2,w_chm_y2,w_chm_z2),logic_w_chm2,"W_chamber2",logic_w,false,0);
        
        // Al Insert#2: tungsten chamber
	G4double w_chm_d_Al2=0.404*mm; // depth
	G4double w_chm_hd_Al2=0.5*w_chm_d_Al2*mm; // half depth
	G4double w_chm_ir_Al2=30.2*mm; // inner radius
	G4double w_chm_or_Al2=33.7*mm; // outer radius
	G4double w_chm_sta_Al2=0.0*deg; // start angle
	G4double w_chm_spa_Al2=360*deg; // span angle
	G4double w_chm_x_Al2=0.0*mm; // x location
	G4double w_chm_y_Al2=0.0*mm; // y location
	G4double w_chm_z_Al2=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_d_Al1+w_chm_d2+w_chm_hd_Al2; // z location
	G4Tubs* solid_w_chm_Al2 = new G4Tubs("W_chamber_Al2",w_chm_ir_Al2,w_chm_or_Al2,w_chm_hd_Al2,w_chm_sta_Al2,w_chm_spa_Al2);
        logic_w_chm_Al2 = new G4LogicalVolume(solid_w_chm_Al2,W,"W_chamber_Al2",0,0,0);
        physi_w_chm_Al2 = new G4PVPlacement(0,G4ThreeVector(w_chm_x_Al2,w_chm_y_Al2,w_chm_z_Al2),logic_w_chm_Al2,"W_chamber_Al2",logic_w,false,0);


        // tungsten chamber Part 3
	G4double w_chm_d3=1.8768*mm; // depth
	G4double w_chm_hd3=0.5*w_chm_d3*mm; // half depth
	G4double w_chm_ir3=30.2*mm; // inner radius
	G4double w_chm_or3=33.7*mm; // outer radius
	G4double w_chm_sta3=0.0*deg; // start angle
	G4double w_chm_spa3=360*deg; // span angle
	G4double w_chm_x3=0.0*mm; // x location
	G4double w_chm_y3=0.0*mm; // y location
	G4double w_chm_z3=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_d_Al1+w_chm_d2+w_chm_d_Al2+w_chm_hd3; // z location
	G4Tubs* solid_w_chm3 = new G4Tubs("W_chamber3",w_chm_ir3,w_chm_or3,w_chm_hd3,w_chm_sta3,w_chm_spa3);
        logic_w_chm3 = new G4LogicalVolume(solid_w_chm3,W,"W_chamber3",0,0,0);
        physi_w_chm3 = new G4PVPlacement(0,G4ThreeVector(w_chm_x3,w_chm_y3,w_chm_z3),logic_w_chm3,"W_chamber3",logic_w,false,0);
        
        // Al Insert #3: tungsten chamber
	G4double w_chm_d_Al3=0.404*mm; // depth
	G4double w_chm_hd_Al3=0.5*w_chm_d_Al3*mm; // half depth
	G4double w_chm_ir_Al3=30.2*mm; // inner radius
	G4double w_chm_or_Al3=33.7*mm; // outer radius
	G4double w_chm_sta_Al3=0.0*deg; // start angle
	G4double w_chm_spa_Al3=360*deg; // span angle
	G4double w_chm_x_Al3=0.0*mm; // x location
	G4double w_chm_y_Al3=0.0*mm; // y location
	G4double w_chm_z_Al3=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_d_Al1+w_chm_d2+w_chm_d_Al2+w_chm_d3+w_chm_hd_Al3; // z location
	G4Tubs* solid_w_chm_Al3 = new G4Tubs("W_chamber_Al3",w_chm_ir_Al3,w_chm_or_Al3,w_chm_hd_Al3,w_chm_sta_Al3,w_chm_spa_Al3);
        logic_w_chm_Al3 = new G4LogicalVolume(solid_w_chm_Al3,Alalloy,"W_chamber_Al3",0,0,0);
        physi_w_chm_Al3 = new G4PVPlacement(0,G4ThreeVector(w_chm_x_Al3,w_chm_y_Al3,w_chm_z_Al3),logic_w_chm_Al3,"W_chamber_Al3",logic_w,false,0);

        // tungsten chamber Part 4
	G4double w_chm_d4=1.8768*mm; // depth
	G4double w_chm_hd4=0.5*w_chm_d4*mm; // half depth
	G4double w_chm_ir4=30.2*mm; // inner radius
	G4double w_chm_or4=33.7*mm; // outer radius
	G4double w_chm_sta4=0.0*deg; // start angle
	G4double w_chm_spa4=360*deg; // span angle
	G4double w_chm_x4=0.0*mm; // x location
	G4double w_chm_y4=0.0*mm; // y location
	G4double w_chm_z4=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_d_Al1+w_chm_d2+w_chm_d_Al2+w_chm_d3+w_chm_d_Al3+w_chm_hd4; // z location
	G4Tubs* solid_w_chm4 = new G4Tubs("W_chamber4",w_chm_ir4,w_chm_or4,w_chm_hd4,w_chm_sta4,w_chm_spa4);
        logic_w_chm4 = new G4LogicalVolume(solid_w_chm4,W,"W_chamber4",0,0,0);
        physi_w_chm4 = new G4PVPlacement(0,G4ThreeVector(w_chm_x4,w_chm_y4,w_chm_z4),logic_w_chm4,"W_chamber4",logic_w,false,0);
        
        // Al Insert #4: tungsten chamber 
	G4double w_chm_d_Al4=0.404*mm; // depth
	G4double w_chm_hd_Al4=0.5*w_chm_d_Al4*mm; // half depth
	G4double w_chm_ir_Al4=30.2*mm; // inner radius
	G4double w_chm_or_Al4=33.7*mm; // outer radius
	G4double w_chm_sta_Al4=0.0*deg; // start angle
	G4double w_chm_spa_Al4=360*deg; // span angle
	G4double w_chm_x_Al4=0.0*mm; // x location
	G4double w_chm_y_Al4=0.0*mm; // y location
	G4double w_chm_z_Al4=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_d_Al1+w_chm_d2+w_chm_d_Al2+w_chm_d3+w_chm_d_Al3+w_chm_d4+w_chm_hd_Al4; // z location
	G4Tubs* solid_w_chm_Al4 = new G4Tubs("W_chamber_Al4",w_chm_ir_Al4,w_chm_or_Al4,w_chm_hd_Al4,w_chm_sta_Al4,w_chm_spa_Al4);
        logic_w_chm_Al4 = new G4LogicalVolume(solid_w_chm_Al4,Alalloy,"W_chamber_Al4",0,0,0);
        physi_w_chm_Al4 = new G4PVPlacement(0,G4ThreeVector(w_chm_x_Al4,w_chm_y_Al4,w_chm_z_Al4),logic_w_chm_Al4,"W_chamber_Al4",logic_w,false,0);
        
        // tungsten chamber Part 5
	G4double w_chm_d5=1.8768*mm; // depth
	G4double w_chm_hd5=0.5*w_chm_d5*mm; // half depth
	G4double w_chm_ir5=30.2*mm; // inner radius
	G4double w_chm_or5=33.7*mm; // outer radius
	G4double w_chm_sta5=0.0*deg; // start angle
	G4double w_chm_spa5=360*deg; // span angle
	G4double w_chm_x5=0.0*mm; // x location
	G4double w_chm_y5=0.0*mm; // y location
	G4double w_chm_z5=frontcoll_d+coll_d+al_ann1_d+w_ann_d1+w_ann_d_Al1+w_ann_d2+w_ann_d_Al2+w_ann_d3+w_ann_d_Al3+w_ann_d4+w_chm_d1+w_chm_d_Al1+w_chm_d2+w_chm_d_Al2+w_chm_d3+w_chm_d_Al3+w_chm_d4+w_chm_d_Al4+w_chm_hd5; // z location
	G4Tubs* solid_w_chm5 = new G4Tubs("W_chamber5",w_chm_ir5,w_chm_or5,w_chm_hd5,w_chm_sta5,w_chm_spa5);
        logic_w_chm5 = new G4LogicalVolume(solid_w_chm5,W,"W_chamber5",0,0,0);
        physi_w_chm5 = new G4PVPlacement(0,G4ThreeVector(w_chm_x5,w_chm_y5,w_chm_z5),logic_w_chm5,"W_chamber5",logic_w,false,0);
*/

    // tungsten end plate
	G4double w_end_d=5.0*mm; // depth
	G4double w_end_hd=0.5*w_end_d*mm; // half depth
	G4double w_end_ir=0.0*mm; // inner radius
	G4double w_end_or=34.75*mm; // outer radius
	G4double w_end_sta=0.0*deg; // start angle
	G4double w_end_spa=360*deg; // span angle
	G4double w_end_x=0.0*mm; // x location
	G4double w_end_y=0.0*mm; // y location
	G4double w_end_z= frontcoll_d + coll_d + al_ann1_d + w_ann_d + w_ann_d_Al1 + w_ann_d1 + w_chm_d+w_end_hd; // z location
	G4Tubs* solid_w_end = new G4Tubs("W_end_plate",w_end_ir,w_end_or,w_end_hd,w_end_sta,w_end_spa);
	//G4VSolid* solid_w_end_subfastener1 = new G4SubtractionSolid("W_end_plate_subfast1",solid_w_end,solid_fastener1);
    //    G4VSolid* solid_w_end_subfastener2 = new G4SubtractionSolid("W_end_plate_subfast2",solid_w_end_subfastener1,solid_fastener2);
    //    G4VSolid* solid_w_end_subfastener3 = new G4SubtractionSolid("W_end_plate_subfast3",solid_w_end_subfastener2,solid_fastener3);
        logic_w_end = new G4LogicalVolume(solid_w_end, W,"W_end_plate",0,0,0);
        physi_w_end = new G4PVPlacement(0,G4ThreeVector(w_end_x,w_end_y,w_end_z),logic_w_end,"W_end_plate",logic_w,false,0);
/*
	// tungsten end plate- Part 1
	G4double w_end_d1=0.572*mm; // depth
	G4double w_end_hd1=0.5*w_end_d1*mm; // half depth
	G4double w_end_ir1=0.0*mm; // inner radius
	G4double w_end_or1=33.7*mm; // outer radius
	G4double w_end_sta1=0.0*deg; // start angle
	G4double w_end_spa1=360*deg; // span angle
	G4double w_end_x1=0.0*mm; // x location
	G4double w_end_y1=0.0*mm; // y location
	G4double w_end_z1=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_hd1; // z location
	G4Tubs* solid_w_end1 = new G4Tubs("W_end_plate1",w_end_ir1,w_end_or1,w_end_hd1,w_end_sta1,w_end_spa1);
	G4VSolid* solid_w_end1_subfastener1 = new G4SubtractionSolid("W_end1_plate_subfast1",solid_w_end1,solid_fastener1);
        G4VSolid* solid_w_end1_subfastener2 = new G4SubtractionSolid("W_end1_plate_subfast2",solid_w_end1_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end1_subfastener3 = new G4SubtractionSolid("W_end1_plate_subfast3",solid_w_end1_subfastener2,solid_fastener3);
        logic_w_end1 = new G4LogicalVolume(solid_w_end1_subfastener3,W,"W_end1_plate_subfast3",0,0,0);
        physi_w_end1 = new G4PVPlacement(0,G4ThreeVector(w_end_x1,w_end_y1,w_end_z1),logic_w_end1,"W_end1_plate_subfast3",logic_w,false,0);
        
        // Al Insert#1: tungsten end plate
	G4double w_end_d_Al1=0.404*mm; // depth
	G4double w_end_hd_Al1=0.5*w_end_d_Al1*mm; // half depth
	G4double w_end_ir_Al1=0.0*mm; // inner radius
	G4double w_end_or_Al1=33.7*mm; // outer radius
	G4double w_end_sta_Al1=0.0*deg; // start angle
	G4double w_end_spa_Al1=360*deg; // span angle
	G4double w_end_x_Al1=0.0*mm; // x location
	G4double w_end_y_Al1=0.0*mm; // y location
	G4double w_end_z_Al1=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_hd_Al1; // z location
	G4Tubs* solid_w_end_Al1 = new G4Tubs("W_end_plate_Al1",w_end_ir_Al1,w_end_or_Al1,w_end_hd_Al1,w_end_sta_Al1,w_end_spa_Al1);
	G4VSolid* solid_w_end_Al1_subfastener1 = new G4SubtractionSolid("W_end_Al1_plate_subfast1",solid_w_end_Al1,solid_fastener1);
        G4VSolid* solid_w_end_Al1_subfastener2 = new G4SubtractionSolid("W_end_Al1_plate_subfast2",solid_w_end_Al1_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end_Al1_subfastener3 = new G4SubtractionSolid("W_end_Al1_plate_subfast3",solid_w_end_Al1_subfastener2,solid_fastener3);
        logic_w_end_Al1 = new G4LogicalVolume(solid_w_end_Al1_subfastener3,Alalloy,"W_end_Al1_plate_subfast3",0,0,0);
        physi_w_end_Al1 = new G4PVPlacement(0,G4ThreeVector(w_end_x_Al1,w_end_y_Al1,w_end_z_Al1),logic_w_end_Al1,"W_end_Al1_plate_subfast3",logic_w,false,0);
       
       	// tungsten end plate- Part 2
	G4double w_end_d2=0.572*mm; // depth
	G4double w_end_hd2=0.5*w_end_d2*mm; // half depth
	G4double w_end_ir2=0.0*mm; // inner radius
	G4double w_end_or2=33.7*mm; // outer radius
	G4double w_end_sta2=0.0*deg; // start angle
	G4double w_end_spa2=360*deg; // span angle
	G4double w_end_x2=0.0*mm; // x location
	G4double w_end_y2=0.0*mm; // y location
	G4double w_end_z2=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_d_Al1+w_end_hd2; // z location
	G4Tubs* solid_w_end2 = new G4Tubs("W_end_plate2",w_end_ir2,w_end_or2,w_end_hd2,w_end_sta2,w_end_spa2);
	G4VSolid* solid_w_end2_subfastener1 = new G4SubtractionSolid("W_end2_plate_subfast1",solid_w_end2,solid_fastener1);
        G4VSolid* solid_w_end2_subfastener2 = new G4SubtractionSolid("W_end2_plate_subfast2",solid_w_end2_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end2_subfastener3 = new G4SubtractionSolid("W_end2_plate_subfast3",solid_w_end2_subfastener2,solid_fastener3);
        logic_w_end2 = new G4LogicalVolume(solid_w_end2_subfastener3,W,"W_end2_plate_subfast3",0,0,0);
        physi_w_end2 = new G4PVPlacement(0,G4ThreeVector(w_end_x2,w_end_y2,w_end_z2),logic_w_end2,"W_end2_plate_subfast3",logic_w,false,0);
        
        // Al Insert#2: tungsten end plate
	G4double w_end_d_Al2=0.404*mm; // depth
	G4double w_end_hd_Al2=0.5*w_end_d_Al2*mm; // half depth
	G4double w_end_ir_Al2=0.0*mm; // inner radius
	G4double w_end_or_Al2=33.7*mm; // outer radius
	G4double w_end_sta_Al2=0.0*deg; // start angle
	G4double w_end_spa_Al2=360*deg; // span angle
	G4double w_end_x_Al2=0.0*mm; // x location
	G4double w_end_y_Al2=0.0*mm; // y location
	G4double w_end_z_Al2=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_d_Al1+w_end_d2+w_end_hd_Al2; // z location
	G4Tubs* solid_w_end_Al2 = new G4Tubs("W_end_plate_Al2",w_end_ir_Al2,w_end_or_Al2,w_end_hd_Al2,w_end_sta_Al2,w_end_spa_Al2);
	G4VSolid* solid_w_end_Al2_subfastener1 = new G4SubtractionSolid("W_end_Al2_plate_subfast1",solid_w_end_Al2,solid_fastener1);
        G4VSolid* solid_w_end_Al2_subfastener2 = new G4SubtractionSolid("W_end_Al2_plate_subfast2",solid_w_end_Al2_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end_Al2_subfastener3 = new G4SubtractionSolid("W_end_Al2_plate_subfast3",solid_w_end_Al2_subfastener2,solid_fastener3);
        logic_w_end_Al2 = new G4LogicalVolume(solid_w_end_Al2_subfastener3,Alalloy,"W_end_Al2_plate_subfast3",0,0,0);
        physi_w_end_Al2 = new G4PVPlacement(0,G4ThreeVector(w_end_x_Al2,w_end_y_Al2,w_end_z_Al2),logic_w_end_Al2,"W_end_Al2_plate_subfast3",logic_w,false,0);

	// tungsten end plate- Part 3
	G4double w_end_d3=0.572*mm; // depth
	G4double w_end_hd3=0.5*w_end_d3*mm; // half depth
	G4double w_end_ir3=0.0*mm; // inner radius
	G4double w_end_or3=33.7*mm; // outer radius
	G4double w_end_sta3=0.0*deg; // start angle
	G4double w_end_spa3=360*deg; // span angle
	G4double w_end_x3=0.0*mm; // x location
	G4double w_end_y3=0.0*mm; // y location
	G4double w_end_z3=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_d_Al1+w_end_d2+w_end_d_Al2+w_end_hd3; // z location
	G4Tubs* solid_w_end3 = new G4Tubs("W_end_plate3",w_end_ir3,w_end_or3,w_end_hd3,w_end_sta3,w_end_spa3);
	G4VSolid* solid_w_end3_subfastener1 = new G4SubtractionSolid("W_end3_plate_subfast1",solid_w_end3,solid_fastener1);
        G4VSolid* solid_w_end3_subfastener2 = new G4SubtractionSolid("W_end3_plate_subfast2",solid_w_end3_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end3_subfastener3 = new G4SubtractionSolid("W_end3_plate_subfast3",solid_w_end3_subfastener2,solid_fastener3);
        logic_w_end3 = new G4LogicalVolume(solid_w_end3_subfastener3,W,"W_end3_plate_subfast3",0,0,0);
        physi_w_end3 = new G4PVPlacement(0,G4ThreeVector(w_end_x3,w_end_y3,w_end_z3),logic_w_end3,"W_end3_plate_subfast3",logic_w,false,0);
        
        // Al Insert#3: tungsten end plate
	G4double w_end_d_Al3=0.404*mm; // depth
	G4double w_end_hd_Al3=0.5*w_end_d_Al3*mm; // half depth
	G4double w_end_ir_Al3=0.0*mm; // inner radius
	G4double w_end_or_Al3=33.7*mm; // outer radius
	G4double w_end_sta_Al3=0.0*deg; // start angle
	G4double w_end_spa_Al3=360*deg; // span angle
	G4double w_end_x_Al3=0.0*mm; // x location
	G4double w_end_y_Al3=0.0*mm; // y location
	G4double w_end_z_Al3=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_d_Al1+w_end_d2+w_end_d_Al2+w_end_d3+w_end_hd_Al3; // z location
	G4Tubs* solid_w_end_Al3 = new G4Tubs("W_end_plate_Al3",w_end_ir_Al3,w_end_or_Al3,w_end_hd_Al3,w_end_sta_Al3,w_end_spa_Al3);
	G4VSolid* solid_w_end_Al3_subfastener1 = new G4SubtractionSolid("W_end_Al3_plate_subfast1",solid_w_end_Al3,solid_fastener1);
        G4VSolid* solid_w_end_Al3_subfastener2 = new G4SubtractionSolid("W_end_Al3_plate_subfast2",solid_w_end_Al3_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end_Al3_subfastener3 = new G4SubtractionSolid("W_end_Al3_plate_subfast3",solid_w_end_Al3_subfastener2,solid_fastener3);
        logic_w_end_Al3 = new G4LogicalVolume(solid_w_end_Al3_subfastener3,Alalloy,"W_end_Al3_plate_subfast3",0,0,0);
        physi_w_end_Al3 = new G4PVPlacement(0,G4ThreeVector(w_end_x_Al3,w_end_y_Al3,w_end_z_Al3),logic_w_end_Al3,"W_end_Al3_plate_subfast3",logic_w,false,0);

	// tungsten end plate- Part 4
	G4double w_end_d4=0.572*mm; // depth
	G4double w_end_hd4=0.5*w_end_d4*mm; // half depth
	G4double w_end_ir4=0.0*mm; // inner radius
	G4double w_end_or4=33.7*mm; // outer radius
	G4double w_end_sta4=0.0*deg; // start angle
	G4double w_end_spa4=360*deg; // span angle
	G4double w_end_x4=0.0*mm; // x location
	G4double w_end_y4=0.0*mm; // y location
	G4double w_end_z4=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_d_Al1+w_end_d2+w_end_d_Al2+w_end_d3+w_end_d_Al3+w_end_hd4; // z location
	G4Tubs* solid_w_end4 = new G4Tubs("W_end_plate4",w_end_ir4,w_end_or4,w_end_hd4,w_end_sta4,w_end_spa4);
	G4VSolid* solid_w_end4_subfastener1 = new G4SubtractionSolid("W_end4_plate_subfast1",solid_w_end4,solid_fastener1);
        G4VSolid* solid_w_end4_subfastener2 = new G4SubtractionSolid("W_end4_plate_subfast2",solid_w_end4_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end4_subfastener3 = new G4SubtractionSolid("W_end4_plate_subfast3",solid_w_end4_subfastener2,solid_fastener3);
        logic_w_end4 = new G4LogicalVolume(solid_w_end4_subfastener3,W,"W_end4_plate_subfast3",0,0,0);
        physi_w_end4 = new G4PVPlacement(0,G4ThreeVector(w_end_x4,w_end_y4,w_end_z4),logic_w_end4,"W_end4_plate_subfast3",logic_w,false,0);
        
        // Al Insert#4: tungsten end plate
	G4double w_end_d_Al4=0.404*mm; // depth
	G4double w_end_hd_Al4=0.5*w_end_d_Al4*mm; // half depth
	G4double w_end_ir_Al4=0.0*mm; // inner radius
	G4double w_end_or_Al4=33.7*mm; // outer radius
	G4double w_end_sta_Al4=0.0*deg; // start angle
	G4double w_end_spa_Al4=360*deg; // span angle
	G4double w_end_x_Al4=0.0*mm; // x location
	G4double w_end_y_Al4=0.0*mm; // y location
	G4double w_end_z_Al4=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_d_Al1+w_end_d2+w_end_d_Al2+w_end_d3+w_end_d_Al3+w_end_d4+w_end_hd_Al4; // z location
	G4Tubs* solid_w_end_Al4 = new G4Tubs("W_end_plate_Al4",w_end_ir_Al4,w_end_or_Al4,w_end_hd_Al4,w_end_sta_Al4,w_end_spa_Al4);
	G4VSolid* solid_w_end_Al4_subfastener1 = new G4SubtractionSolid("W_end_Al4_plate_subfast1",solid_w_end_Al4,solid_fastener1);
        G4VSolid* solid_w_end_Al4_subfastener2 = new G4SubtractionSolid("W_end_Al4_plate_subfast2",solid_w_end_Al4_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end_Al4_subfastener3 = new G4SubtractionSolid("W_end_Al4_plate_subfast3",solid_w_end_Al4_subfastener2,solid_fastener3);
        logic_w_end_Al4 = new G4LogicalVolume(solid_w_end_Al4_subfastener3,Alalloy,"W_end_Al4_plate_subfast3",0,0,0);
        physi_w_end_Al4 = new G4PVPlacement(0,G4ThreeVector(w_end_x_Al4,w_end_y_Al4,w_end_z_Al4),logic_w_end_Al4,"W_end_Al4_plate_subfast3",logic_w,false,0);

	// tungsten end plate- Part 5
	G4double w_end_d5=0.572*mm; // depth
	G4double w_end_hd5=0.5*w_end_d5*mm; // half depth
	G4double w_end_ir5=0.0*mm; // inner radius
	G4double w_end_or5=33.7*mm; // outer radius
	G4double w_end_sta5=0.0*deg; // start angle
	G4double w_end_spa5=360*deg; // span angle
	G4double w_end_x5=0.0*mm; // x location
	G4double w_end_y5=0.0*mm; // y location
	G4double w_end_z5=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d1+w_end_d_Al1+w_end_d2+w_end_d_Al2+w_end_d3+w_end_d_Al3+w_end_d4+w_end_d_Al4+w_end_hd5; // z location
	G4Tubs* solid_w_end5 = new G4Tubs("W_end_plate5",w_end_ir5,w_end_or5,w_end_hd5,w_end_sta5,w_end_spa5);
	G4VSolid* solid_w_end5_subfastener1 = new G4SubtractionSolid("W_end5_plate_subfast1",solid_w_end5,solid_fastener1);
        G4VSolid* solid_w_end5_subfastener2 = new G4SubtractionSolid("W_end5_plate_subfast2",solid_w_end5_subfastener1,solid_fastener2);
        G4VSolid* solid_w_end5_subfastener3 = new G4SubtractionSolid("W_end5_plate_subfast3",solid_w_end5_subfastener2,solid_fastener3);
        logic_w_end5 = new G4LogicalVolume(solid_w_end5_subfastener3,W,"W_end5_plate_subfast3",0,0,0);
        physi_w_end5 = new G4PVPlacement(0,G4ThreeVector(w_end_x5,w_end_y5,w_end_z5),logic_w_end5,"W_end5_plate_subfast3",logic_w,false,0);
        
        // tangsten end enhancing plate
	//G4double w_endenh_d=2.0*mm; // depth
	//G4double w_endenh_hd=0.5*w_endenh_d*mm; // half depth
	//G4double w_endenh_ir=0.0*mm; // inner radius
	//G4double w_endenh_or=30.2*mm; // outer radius
	//G4double w_endenh_sta=0.0*deg; // start angle
	//G4double w_endenh_spa=360*deg; // span angle
	//G4double w_endenh_x=0.0*mm; // x location
	//G4double w_endenh_y=0.0*mm; // y location
	//G4double w_endenh_z=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d+w_endenh_hd; // z location
	//G4Tubs* solid_w_endenh = new G4Tubs("W_endenh_plate",w_endenh_ir,w_endenh_or,w_endenh_hd,w_endenh_sta,w_endenh_spa);
 //       G4VSolid* solid_w_endenh_subfastener1 = new G4SubtractionSolid("W_endenh_plate_subfast1",solid_w_endenh,solid_fastener1);
 //       G4VSolid* solid_w_endenh_subfastener2 = new G4SubtractionSolid("W_endenh_plate_subfast2",solid_w_endenh_subfastener1,solid_fastener2);
 //       G4VSolid* solid_w_endenh_subfastener3 = new G4SubtractionSolid("W_endenh_plate_subfast3",solid_w_endenh_subfastener2,solid_fastener3);
 //       logic_w_endenh = new G4LogicalVolume(solid_w_endenh_subfastener3,W,"W_endenh_plate_subfast3",0,0,0);
 //       physi_w_endenh = new G4PVPlacement(0,G4ThreeVector(w_endenh_x,w_endenh_y,w_endenh_z),logic_w_endenh,"W_endenh_plate_subfast3",logic_w,false,0);
	
	// tangsten end enhancing plate - Part 1
	G4double w_endenh_d1=2.0*mm; // depth
	G4double w_endenh_hd1=0.5*w_endenh_d1*mm; // half depth
	G4double w_endenh_ir1=0.0*mm; // inner radius
	G4double w_endenh_or1=30.2*mm; // outer radius
	G4double w_endenh_sta1=0.0*deg; // start angle
	G4double w_endenh_spa1=360*deg; // span angle
	G4double w_endenh_x1=0.0*mm; // x location
	G4double w_endenh_y1=0.0*mm; // y location
	G4double w_endenh_z1=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d+w_endenh_hd1; // z location
	G4Tubs* solid_w_endenh1 = new G4Tubs("W_endenh1_plate",w_endenh_ir1,w_endenh_or1,w_endenh_hd1,w_endenh_sta1,w_endenh_spa1);
        G4VSolid* solid_w_endenh1_subfastener1 = new G4SubtractionSolid("W_endenh1_plate_subfast1",solid_w_endenh1,solid_fastener1);
        G4VSolid* solid_w_endenh1_subfastener2 = new G4SubtractionSolid("W_endenh1_plate_subfast2",solid_w_endenh1_subfastener1,solid_fastener2);
        G4VSolid* solid_w_endenh1_subfastener3 = new G4SubtractionSolid("W_endenh1_plate_subfast3",solid_w_endenh1_subfastener2,solid_fastener3);
        logic_w_endenh1 = new G4LogicalVolume(solid_w_endenh1_subfastener3,W,"W_endenh1_plate_subfast3",0,0,0);
        physi_w_endenh1 = new G4PVPlacement(0,G4ThreeVector(w_endenh_x1,w_endenh_y1,w_endenh_z1),logic_w_endenh1,"W_endenh1_plate_subfast3",logic_w,false,0);

         // Al Insert #1: tangsten end enhancing plate
	G4double w_endenh_d_Al1=2.0*mm; // depth
	G4double w_endenh_hd_Al1=0.5*w_endenh_d_Al1*mm; // half depth
	G4double w_endenh_ir_Al1=0.0*mm; // inner radius
	G4double w_endenh_or_Al1=30.2*mm; // outer radius
	G4double w_endenh_sta_Al1=0.0*deg; // start angle
	G4double w_endenh_spa_Al1=360*deg; // span angle
	G4double w_endenh_x_Al1=0.0*mm; // x location
	G4double w_endenh_y_Al1=0.0*mm; // y location
	G4double w_endenh_z_Al1=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d+w_endenh_d1+w_endenh_hd_Al1; // z location
	G4Tubs* solid_w_endenh_Al1 = new G4Tubs("W_endenh_Al1_plate",w_endenh_ir_Al1,w_endenh_or_Al1,w_endenh_hd_Al1,w_endenh_sta_Al1,w_endenh_spa_Al1);
        G4VSolid* solid_w_endenh_Al1_subfastener1 = new G4SubtractionSolid("W_endenh_Al1_plate_subfast1",solid_w_endenh_Al1,solid_fastener1);
        G4VSolid* solid_w_endenh_Al1_subfastener2 = new G4SubtractionSolid("W_endenh_Al1_plate_subfast2",solid_w_endenh_Al1_subfastener1,solid_fastener2);
        G4VSolid* solid_w_endenh_Al1_subfastener3 = new G4SubtractionSolid("W_endenh_Al1_plate_subfast3",solid_w_endenh_Al1_subfastener2,solid_fastener3);
        logic_w_endenh_Al1 = new G4LogicalVolume(solid_w_endenh_Al1_subfastener3,Alalloy,"W_endenh_Al1_plate_subfast3",0,0,0);
        physi_w_endenh_Al1 = new G4PVPlacement(0,G4ThreeVector(w_endenh_x_Al1,w_endenh_y_Al1,w_endenh_z_Al1),logic_w_endenh_Al1,"W_endenh_Al1_plate_subfast3",logic_w,false,0);
       
       // tangsten end enhancing plate - Part 2
	G4double w_endenh_d2=2.0*mm; // depth
	G4double w_endenh_hd2=0.5*w_endenh_d2*mm; // half depth
	G4double w_endenh_ir2=0.0*mm; // inner radius
	G4double w_endenh_or2=30.2*mm; // outer radius
	G4double w_endenh_sta2=0.0*deg; // start angle
	G4double w_endenh_spa2=360*deg; // span angle
	G4double w_endenh_x2=0.0*mm; // x location
	G4double w_endenh_y2=0.0*mm; // y location
	G4double w_endenh_z2=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d+w_endenh_d1+w_endenh_d_Al1+w_endenh_hd2; // z location
	G4Tubs* solid_w_endenh2 = new G4Tubs("W_endenh2_plate",w_endenh_ir2,w_endenh_or2,w_endenh_hd2,w_endenh_sta2,w_endenh_spa2);
        G4VSolid* solid_w_endenh2_subfastener1 = new G4SubtractionSolid("W_endenh2_plate_subfast1",solid_w_endenh2,solid_fastener1);
        G4VSolid* solid_w_endenh2_subfastener2 = new G4SubtractionSolid("W_endenh2_plate_subfast2",solid_w_endenh2_subfastener1,solid_fastener2);
        G4VSolid* solid_w_endenh2_subfastener3 = new G4SubtractionSolid("W_endenh2_plate_subfast3",solid_w_endenh2_subfastener2,solid_fastener3);
        logic_w_endenh2 = new G4LogicalVolume(solid_w_endenh2_subfastener3,W,"W_endenh2_plate_subfast3",0,0,0);
        physi_w_endenh2 = new G4PVPlacement(0,G4ThreeVector(w_endenh_x2,w_endenh_y2,w_endenh_z2),logic_w_endenh2,"W_endenh2_plate_subfast3",logic_w,false,0);

         // Al Insert #2: tangsten end enhancing plate
	G4double w_endenh_d_Al2=2.0*mm; // depth
	G4double w_endenh_hd_Al2=0.5*w_endenh_d_Al2*mm; // half depth
	G4double w_endenh_ir_Al2=0.0*mm; // inner radius
	G4double w_endenh_or_Al2=30.2*mm; // outer radius
	G4double w_endenh_sta_Al2=0.0*deg; // start angle
	G4double w_endenh_spa_Al2=360*deg; // span angle
	G4double w_endenh_x_Al2=0.0*mm; // x location
	G4double w_endenh_y_Al2=0.0*mm; // y location
	G4double w_endenh_z_Al2=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d+w_endenh_d1+w_endenh_d_Al1+w_endenh_d2+w_endenh_hd_Al2; // z location
	G4Tubs* solid_w_endenh_Al2 = new G4Tubs("W_endenh_Al2_plate",w_endenh_ir_Al2,w_endenh_or_Al2,w_endenh_hd_Al2,w_endenh_sta_Al2,w_endenh_spa_Al2);
        G4VSolid* solid_w_endenh_Al2_subfastener1 = new G4SubtractionSolid("W_endenh_Al2_plate_subfast1",solid_w_endenh_Al2,solid_fastener1);
        G4VSolid* solid_w_endenh_Al2_subfastener2 = new G4SubtractionSolid("W_endenh_Al2_plate_subfast2",solid_w_endenh_Al2_subfastener1,solid_fastener2);
        G4VSolid* solid_w_endenh_Al2_subfastener3 = new G4SubtractionSolid("W_endenh_Al2_plate_subfast3",solid_w_endenh_Al2_subfastener2,solid_fastener3);
        logic_w_endenh_Al2 = new G4LogicalVolume(solid_w_endenh_Al2_subfastener3,Alalloy,"W_endenh_Al2_plate_subfast3",0,0,0);
        physi_w_endenh_Al2 = new G4PVPlacement(0,G4ThreeVector(w_endenh_x_Al2,w_endenh_y_Al2,w_endenh_z_Al2),logic_w_endenh_Al2,"W_endenh_Al2_plate_subfast3",logic_w,false,0);
        
        // tangsten end enhancing plate - Part 3
	G4double w_endenh_d3=2.0*mm; // depth
	G4double w_endenh_hd3=0.5*w_endenh_d3*mm; // half depth
	G4double w_endenh_ir3=0.0*mm; // inner radius
	G4double w_endenh_or3=30.2*mm; // outer radius
	G4double w_endenh_sta3=0.0*deg; // start angle
	G4double w_endenh_spa3=360*deg; // span angle
	G4double w_endenh_x3=0.0*mm; // x location
	G4double w_endenh_y3=0.0*mm; // y location
	G4double w_endenh_z3=frontcoll_d+coll_d+al_ann1_d+w_ann_d+w_chm_d+w_end_d+w_endenh_d1+w_endenh_d_Al1+w_endenh_d2+w_endenh_d_Al2+w_endenh_hd3; // z location
	G4Tubs* solid_w_endenh3 = new G4Tubs("W_endenh3_plate",w_endenh_ir3,w_endenh_or3,w_endenh_hd3,w_endenh_sta3,w_endenh_spa3);
        G4VSolid* solid_w_endenh3_subfastener1 = new G4SubtractionSolid("W_endenh3_plate_subfast1",solid_w_endenh3,solid_fastener1);
        G4VSolid* solid_w_endenh3_subfastener2 = new G4SubtractionSolid("W_endenh3_plate_subfast2",solid_w_endenh3_subfastener1,solid_fastener2);
        G4VSolid* solid_w_endenh3_subfastener3 = new G4SubtractionSolid("W_endenh3_plate_subfast3",solid_w_endenh3_subfastener2,solid_fastener3);
        logic_w_endenh3 = new G4LogicalVolume(solid_w_endenh3_subfastener3,W,"W_endenh3_plate_subfast3",0,0,0);
        physi_w_endenh3 = new G4PVPlacement(0,G4ThreeVector(w_endenh_x3,w_endenh_y3,w_endenh_z3),logic_w_endenh3,"W_endenh3_plate_subfast3",logic_w,false,0);
	*/

	// detectors (d)
	// detector 1 (d1)
	G4double d1_d = 1.5*mm; // depth
	G4double d1_hd = 0.5*d1_d*mm; // half depth
	G4double d1_ir = 0.0*mm; // inner radius
	G4double d1_or = 10.0*mm; // outer radius
	G4double d1_sta = 0.0*deg; // start angle
	G4double d1_spa = 360*deg; // span angle
	G4double d1_x=0.0*mm; // x location
	G4double d1_y=0.0*mm; // y location
	G4double d1_z= w_ann_z1+ w_ann_hd1+0.5*mm + d1_hd; // z location
	G4Tubs* solid_d1 = new G4Tubs("detector_1",d1_ir,d1_or,d1_hd,d1_sta,d1_spa);
	logic_d1 = new G4LogicalVolume(solid_d1,Si,"detector_1",0,0,0);
	physi_d1 = new G4PVPlacement(0,G4ThreeVector(d1_x,d1_y,d1_z),logic_d1,"detector_1",logic_w,false,0);

	//Source Spherical Cap
	G4double pRmin = 85*mm;// *mm;
	G4double pRmax = pRmin + 0.2 * mm;
	G4double pSPhi = 0 * deg;
	G4double pDPhi = 360 * deg;
	G4double pSTheta = 165* deg;
	G4double pDTheta = 180 * deg - pSTheta;
	G4double S1_x = 0.0 * mm; // x location
	G4double S1_y = 0.0 * mm; // y location
	G4double S1_z = d1_z; // z location-centered on first detector  d1_z
	const G4String& pName = "source";

	G4Sphere* S1 = new G4Sphere(pName, pRmin, pRmax, pSPhi, pDPhi, pSTheta, pDTheta);
	logic_S1 = new G4LogicalVolume(S1, Vacuum, "source", 0, 0, 0);
	physi_S1 = new G4PVPlacement(0, G4ThreeVector(S1_x, S1_y, S1_z), logic_S1, "source", logic_w, false, 0);

        /* // detector 1 outer (d1)
        G4double d1outer_d = 1.5*mm; // depth
        G4double d1outer_hd = 0.5*d1outer_d*mm; // half depth
        G4double d1outer_ir = 10.0*mm; // inner radius
        G4double d1outer_or = 20.0*mm; // outer radius
        G4double d1outer_sta = 0.0*deg; // start angle
        G4double d1outer_spa = 360*deg; // span angle
        G4double d1outer_x = 0.0*mm; // x location
        G4double d1outer_y = 0.0*mm; // y location
        G4double d1outer_z=frontcoll_d+coll_d+al_ann1_d+w_ann_d+d1_hd+1.0*mm; // z location
        G4Tubs* solid_d1outer = new G4Tubs("detector_1Outer",d1outer_ir,d1outer_or,d1outer_hd,d1outer_sta,d1outer_spa);
        logic_d1outer = new G4LogicalVolume(solid_d1outer,Si,"detector_1Outer",0,0,0);
        physi_d1outer = new G4PVPlacement(0,G4ThreeVector(d1outer_x,d1outer_y,d1outer_z),logic_d1outer,"detector_1Outer",logic_w,false,0); */
        
   	// detector 2 (d2)
        G4double d2_d = 1.5*mm; // depth
        G4double d2_hd = 0.5*d2_d*mm; // half depth
        G4double d2_ir = 0.0*mm; // inner radius
        G4double d2_or = 20.0*mm; // outer radius
        G4double d2_sta = 0.0*deg; // start angle
        G4double d2_spa = 360*deg; // span angle
        G4double d2_x=0.0*mm; // x location
        G4double d2_y=0.0*mm; // y location
        G4double d2_z= d1_z+ d1_hd +1.080*mm+d2_hd; // z location
        G4Tubs* solid_d2 = new G4Tubs("detector_2",d2_ir,d2_or,d2_hd,d2_sta,d2_spa);
        logic_d2 = new G4LogicalVolume(solid_d2,Si,"detector_2",0,0,0);
        physi_d2 = new G4PVPlacement(0,G4ThreeVector(d2_x,d2_y,d2_z),logic_d2,"detector_2",logic_w,false,0);
    
        /*/ detector 2 outer(d2)
        G4double d2outer_d = 1.5*mm; // depth
        G4double d2outer_hd = 0.5*d2outer_d*mm; // half depth
        G4double d2outer_ir = 10.0*mm; // inner radius
        G4double d2outer_or = 20.0*mm; // outer radius
        G4double d2outer_sta = 0.0*deg; // start angle
        G4double d2outer_spa = 360*deg; // span angle
        G4double d2outer_x = 0.0*mm; // x location
        G4double d2outer_y = 0.0*mm; // y location
        G4double d2outer_z = frontcoll_d+coll_d+al_ann1_d+w_ann_d+d1_d+d2outer_hd+2.0*mm; // z location
        G4Tubs* solid_d2outer = new G4Tubs("detector_2Outer",d2outer_ir,d2outer_or,d2outer_hd,d2outer_sta,d2outer_spa);
        logic_d2outer = new G4LogicalVolume(solid_d2outer,Si,"detector_2Outer",0,0,0);
        physi_d2outer = new G4PVPlacement(0,G4ThreeVector(d2outer_x,d2outer_y,d2outer_z),logic_d2outer,"detector_2Outer",logic_w,false,0);*/
        
	// detector 3 (d3)
	G4double d3_d = 1.5*mm; // depth
	G4double d3_hd = 0.5*d3_d*mm; // half depth
	G4double d3_ir = 0.0*mm; // inner radius
	G4double d3_or = 20.0*mm; // outer radius
	G4double d3_sta = 0.0*deg; // start angle
	G4double d3_spa = 360*deg; // span angle
	G4double d3_x=0.0*mm; // x location
	G4double d3_y=0.0*mm; // y location
	G4double d3_z= d2_z + d2_hd + 1.00 * mm + d3_hd; // z location
	G4Tubs* solid_d3 = new G4Tubs("detector_3",d3_ir,d3_or,d3_hd,d3_sta,d3_spa);
	logic_d3 = new G4LogicalVolume(solid_d3,Si,"detector_3",0,0,0);
	physi_d3 = new G4PVPlacement(0,G4ThreeVector(d3_x,d3_y,d3_z),logic_d3,"detector_3",logic_w,false,0);

        /*// detector 3 outer (d3)
        G4double d3outer_d = 1.5*mm; // depth
        G4double d3outer_hd = 0.5*d3outer_d*mm; // half depth
        G4double d3outer_ir = 10.0*mm; // inner radius
        G4double d3outer_or = 20.0*mm; // outer radius
        G4double d3outer_sta = 0.0*deg; // start angle
        G4double d3outer_spa = 360*deg; // span angle
        G4double d3outer_x = 0.0*mm; // x location
        G4double d3outer_y = 0.0*mm; // y location
        G4double d3outer_z = frontcoll_d+coll_d+al_ann1_d+w_ann_d+d1_d+d2_d+d3_hd+3.0*mm; // z location
        G4Tubs* solid_d3outer = new G4Tubs("detector_3Outer",d3outer_ir,d3outer_or,d3outer_hd,d3outer_sta,d3outer_spa);
        logic_d3outer = new G4LogicalVolume(solid_d3outer,Si,"detector_3Outer",0,0,0);
        physi_d3outer = new G4PVPlacement(0,G4ThreeVector(d3outer_x,d3outer_y,d3outer_z),logic_d3outer,"detector_3Outer",logic_w,false,0);*/
    
	// detector 4 (d4)
	G4double d4_d = 1.5*mm; // depth
	G4double d4_hd = 0.5*d4_d*mm; // half depth
	G4double d4_ir = 0.0*mm; // inner radius
	G4double d4_or = 20.0*mm; // outer radius
	G4double d4_sta = 0.0*deg; // start angle
	G4double d4_spa = 360*deg; // span angle
	G4double d4_x=0.0*mm; // x location
	G4double d4_y=0.0*mm; // y location
	G4double d4_z= d3_z + d3_hd + 1.080 * mm + d4_hd; // z location
	G4Tubs* solid_d4 = new G4Tubs("detector_4",d4_ir,d4_or,d4_hd,d4_sta,d4_spa);
	logic_d4 = new G4LogicalVolume(solid_d4,Si,"detector_4",0,0,0);
	physi_d4 = new G4PVPlacement(0,G4ThreeVector(d4_x,d4_y,d4_z),logic_d4,"detector_4",logic_w,false,0);

        /* detector 4 outer(d4)
        G4double d4outer_d = 1.5*mm; // depth
        G4double d4outer_hd = 0.5*d4outer_d*mm; // half depth
        G4double d4outer_ir = 10.0*mm; // inner radius
        G4double d4outer_or = 20.0*mm; // outer radius
        G4double d4outer_sta = 0.0*deg; // start angle
        G4double d4outer_spa = 360*deg; // span angle
        G4double d4outer_x = 0.0*mm; // x location
        G4double d4outer_y = 0.0*mm; // y location
        G4double d4outer_z = frontcoll_d+coll_d+al_ann1_d+w_ann_d+d1_d+d2_d+d3_d+d4_hd+4.0*mm; // z location
        G4Tubs* solid_d4outer = new G4Tubs("detector_4Outer",d4outer_ir,d4outer_or,d4outer_hd,d4outer_sta,d4outer_spa);
        logic_d4outer = new G4LogicalVolume(solid_d4outer,Si,"detector_4Outer",0,0,0);
        physi_d4outer = new G4PVPlacement(0,G4ThreeVector(d4outer_x,d4outer_y,d4outer_z),logic_d4outer,"detector_4Outer",logic_w,false,0);*/


		// detector 5 (d5)
		G4double d5_d = 1.5 * mm; // depth
		G4double d5_hd = 0.5 * d5_d * mm; // half depth
		G4double d5_ir = 0.0 * mm; // inner radius
		G4double d5_or = 20.0 * mm; // outer radius
		G4double d5_sta = 0.0 * deg; // start angle
		G4double d5_spa = 360 * deg; // span angle
		G4double d5_x = 0.0 * mm; // x location
		G4double d5_y = 0.0 * mm; // y location
		G4double d5_z = d4_z + d4_hd + 1.00 * mm + d5_hd; // z location
		G4Tubs* solid_d5 = new G4Tubs("detector_5", d5_ir, d5_or, d5_hd, d5_sta, d5_spa);
		logic_d5 = new G4LogicalVolume(solid_d5, Si, "detector_5", 0, 0, 0);
		physi_d5 = new G4PVPlacement(0, G4ThreeVector(d5_x, d5_y, d5_z), logic_d5, "detector_5", logic_w, false, 0);

		/*/ detector 5 outer (d5)
		G4double d5outer_d = 1.5 * mm; // depth
		G4double d5outer_hd = 0.5 * d5outer_d * mm; // half depth
		G4double d5outer_ir = 10.0 * mm; // inner radius
		G4double d5outer_or = 20.0 * mm; // outer radius
		G4double d5outer_sta = 0.0 * deg; // start angle
		G4double d5outer_spa = 360 * deg; // span angle
		G4double d5outer_x = 0.0 * mm; // x location
		G4double d5outer_y = 0.0 * mm; // y location
		G4double d5outer_z = frontcoll_d + coll_d + al_ann1_d + w_ann_d + d1_d + d2_d + d3_d + d4_d + d5_hd + 5.0 * mm; // z location
		G4Tubs* solid_d5outer = new G4Tubs("detector_5Outer", d5outer_ir, d5outer_or, d5outer_hd, d5outer_sta, d5outer_spa);
		logic_d5outer = new G4LogicalVolume(solid_d5outer, Si, "detector_5Outer", 0, 0, 0);
		physi_d5outer = new G4PVPlacement(0, G4ThreeVector(d5outer_x, d5outer_y, d5outer_z), logic_d5outer, "detector_5Outer", logic_w, false, 0);*/

		// detector 6 (d6)
		G4double d6_d = 1.5 * mm; // depth
		G4double d6_hd = 0.5 * d6_d * mm; // half depth
		G4double d6_ir = 0.0 * mm; // inner radius
		G4double d6_or = 20.0 * mm; // outer radius
		G4double d6_sta = 0.0 * deg; // start angle
		G4double d6_spa = 360 * deg; // span angle
		G4double d6_x = 0.0 * mm; // x location
		G4double d6_y = 0.0 * mm; // y location
		G4double d6_z = d5_z + d5_hd + 1.080 * mm + d6_hd; // z location
		G4Tubs* solid_d6 = new G4Tubs("detector_6", d6_ir, d6_or, d6_hd, d6_sta, d6_spa);
		logic_d6 = new G4LogicalVolume(solid_d6, Si, "detector_6", 0, 0, 0);
		physi_d6 = new G4PVPlacement(0, G4ThreeVector(d6_x, d6_y, d6_z), logic_d6, "detector_6", logic_w, false, 0);

		/*/ detector 6 outer (d6)
		G4double d6outer_d = 1.5 * mm; // depth
		G4double d6outer_hd = 0.5 * d6outer_d * mm; // half depth
		G4double d6outer_ir = 10.0 * mm; // inner radius
		G4double d6outer_or = 20.0 * mm; // outer radius
		G4double d6outer_sta = 0.0 * deg; // start angle
		G4double d6outer_spa = 360 * deg; // span angle
		G4double d6outer_x = 0.0 * mm; // x location
		G4double d6outer_y = 0.0 * mm; // y location
		G4double d6outer_z = frontcoll_d + coll_d + al_ann1_d + w_ann_d + d1_d + d2_d + d3_d + d4_d + d5_d + d6_hd + 6.0 * mm; // z location
		G4Tubs* solid_d6outer = new G4Tubs("detector_6Outer", d6outer_ir, d6outer_or, d6outer_hd, d6outer_sta, d6outer_spa);
		logic_d6outer = new G4LogicalVolume(solid_d6outer, Si, "detector_6Outer", 0, 0, 0);
		physi_d6outer = new G4PVPlacement(0, G4ThreeVector(d6outer_x, d6outer_y, d6outer_z), logic_d6outer, "detector_6Outer", logic_w, false, 0);*/


		// detector 7 (d7)
		G4double d7_d = 1.5 * mm; // depth
		G4double d7_hd = 0.5 * d7_d * mm; // half depth
		G4double d7_ir = 0.0 * mm; // inner radius
		G4double d7_or = 20.0 * mm; // outer radius
		G4double d7_sta = 0.0 * deg; // start angle
		G4double d7_spa = 360 * deg; // span angle
		G4double d7_x = 0.0 * mm; // x location
		G4double d7_y = 0.0 * mm; // y location
		G4double d7_z = d6_z + d6_hd + 1.00 * mm + d7_hd; // z location
		G4Tubs* solid_d7 = new G4Tubs("detector_7", d7_ir, d7_or, d7_hd, d7_sta, d7_spa);
		logic_d7 = new G4LogicalVolume(solid_d7, Si, "detector_7", 0, 0, 0);
		physi_d7 = new G4PVPlacement(0, G4ThreeVector(d7_x, d7_y, d7_z), logic_d7, "detector_7", logic_w, false, 0);

		/*/ detector 7 outer (d7)
		G4double d7outer_d = 1.5 * mm; // depth
		G4double d7outer_hd = 0.5 * d7outer_d * mm; // half depth
		G4double d7outer_ir = 10.0 * mm; // inner radius
		G4double d7outer_or = 20.0 * mm; // outer radius
		G4double d7outer_sta = 0.0 * deg; // start angle
		G4double d7outer_spa = 360 * deg; // span angle
		G4double d7outer_x = 0.0 * mm; // x location
		G4double d7outer_y = 0.0 * mm; // y location
		G4double d7outer_z = frontcoll_d + coll_d + al_ann1_d + w_ann_d + d1_d + d2_d + d3_d +d4_d + d5_d + d6_d + d7_hd + 7.0 * mm; // z location
		G4Tubs* solid_d7outer = new G4Tubs("detector_7Outer", d7outer_ir, d7outer_or, d7outer_hd, d7outer_sta, d7outer_spa);
		logic_d7outer = new G4LogicalVolume(solid_d7outer, Si, "detector_7Outer", 0, 0, 0);
		physi_d7outer = new G4PVPlacement(0, G4ThreeVector(d7outer_x, d7outer_y, d7outer_z), logic_d7outer, "detector_7Outer", logic_w, false, 0);*/


		// detector 8 (d8)
		G4double d8_d = 1.5 * mm; // depth
		G4double d8_hd = 0.5 * d8_d * mm; // half depth
		G4double d8_ir = 0.0 * mm; // inner radius
		G4double d8_or = 20.0 * mm; // outer radius
		G4double d8_sta = 0.0 * deg; // start angle
		G4double d8_spa = 360 * deg; // span angle
		G4double d8_x = 0.0 * mm; // x location
		G4double d8_y = 0.0 * mm; // y location
		G4double d8_z = d7_z + d7_hd + 1.080 * mm + d8_hd; // z location
		G4Tubs* solid_d8 = new G4Tubs("detector_8", d8_ir, d8_or, d8_hd, d8_sta, d8_spa);
		logic_d8 = new G4LogicalVolume(solid_d8, Si, "detector_8", 0, 0, 0);
		physi_d8 = new G4PVPlacement(0, G4ThreeVector(d8_x, d8_y, d8_z), logic_d8, "detector_8", logic_w, false, 0);

		/*/ detector 8 outer (d8)
		G4double d8outer_d = 1.5 * mm; // depth
		G4double d8outer_hd = 0.5 * d8outer_d * mm; // half depth
		G4double d8outer_ir = 10.0 * mm; // inner radius
		G4double d8outer_or = 20.0 * mm; // outer radius
		G4double d8outer_sta = 0.0 * deg; // start angle
		G4double d8outer_spa = 360 * deg; // span angle
		G4double d8outer_x = 0.0 * mm; // x location
		G4double d8outer_y = 0.0 * mm; // y location
		G4double d8outer_z = frontcoll_d + coll_d + al_ann1_d + w_ann_d + d1_d + d2_d + d3_d + d4_d + d5_d + d6_d + d7_d + d8_hd + 8.0 * mm; // z location
		G4Tubs* solid_d8outer = new G4Tubs("detector_8Outer", d8outer_ir, d8outer_or, d8outer_hd, d8outer_sta, d8outer_spa);
		logic_d8outer = new G4LogicalVolume(solid_d8outer, Si, "detector_8Outer", 0, 0, 0);
		physi_d8outer = new G4PVPlacement(0, G4ThreeVector(d8outer_x, d8outer_y, d8outer_z), logic_d8outer, "detector_8Outer", logic_w, false, 0);*/


		// detector 9 (d9)
		G4double d9_d = 1.5 * mm; // depth
		G4double d9_hd = 0.5 * d9_d * mm; // half depth
		G4double d9_ir = 0.0 * mm; // inner radius
		G4double d9_or = 20.0 * mm; // outer radius
		G4double d9_sta = 0.0 * deg; // start angle
		G4double d9_spa = 360 * deg; // span angle
		G4double d9_x = 0.0 * mm; // x location
		G4double d9_y = 0.0 * mm; // y location
		G4double d9_z = d8_z + d8_hd + 1.00 * mm + d9_hd; // z location
		G4Tubs* solid_d9 = new G4Tubs("detector_9", d9_ir, d9_or, d9_hd, d9_sta, d9_spa);
		logic_d9 = new G4LogicalVolume(solid_d9, Si, "detector_9", 0, 0, 0);
		physi_d9 = new G4PVPlacement(0, G4ThreeVector(d9_x, d9_y, d9_z), logic_d9, "detector_9", logic_w, false, 0);

		// Al Insert: Shims placed behind Detector Stack
		G4double w_end_d_Al1 = 2 * mm; // depth
		G4double w_end_hd_Al1 = 0.5 * w_end_d_Al1 * mm; // half depth
		G4double w_end_ir_Al1 = 0.0 * mm; // inner radius
		G4double w_end_or_Al1 = 29.3 * mm; // outer radius
		G4double w_end_sta_Al1 = 0.0 * deg; // start angle
		G4double w_end_spa_Al1 = 360 * deg; // span angle
		G4double w_end_x_Al1 = 0.0 * mm; // x location
		G4double w_end_y_Al1 = 0.0 * mm; // y location
		G4double w_end_z_Al1 = d9_z + d9_hd +0.54*mm+ w_end_hd_Al1; // z location
		G4Tubs* solid_w_end_Al1 = new G4Tubs("W_end_plate_Al1", w_end_ir_Al1, w_end_or_Al1, w_end_hd_Al1, w_end_sta_Al1, w_end_spa_Al1);
		//G4VSolid* solid_w_end_Al1_subfastener1 = new G4SubtractionSolid("W_end_Al1_plate_subfast1", solid_w_end_Al1, solid_fastener1);
		//G4VSolid* solid_w_end_Al1_subfastener2 = new G4SubtractionSolid("W_end_Al1_plate_subfast2", solid_w_end_Al1_subfastener1, solid_fastener2);
		//G4VSolid* solid_w_end_Al1_subfastener3 = new G4SubtractionSolid("W_end_Al1_plate_subfast3", solid_w_end_Al1_subfastener2, solid_fastener3);
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
		//G4VSolid* solid_w_end_Al1_subfastener1 = new G4SubtractionSolid("W_end_Al1_plate_subfast1", solid_w_end_Al1, solid_fastener1);
		//G4VSolid* solid_w_end_Al1_subfastener2 = new G4SubtractionSolid("W_end_Al1_plate_subfast2", solid_w_end_Al1_subfastener1, solid_fastener2);
		//G4VSolid* solid_w_end_Al1_subfastener3 = new G4SubtractionSolid("W_end_Al1_plate_subfast3", solid_w_end_Al1_subfastener2, solid_fastener3);
		logic_w_end_Al2 = new G4LogicalVolume(solid_w_end_Al2, Alalloy, "W_end_Al2", 0, 0, 0);
		physi_w_end_Al2 = new G4PVPlacement(0, G4ThreeVector(w_end_x_Al2, w_end_y_Al2, w_end_z_Al2), logic_w_end_Al2, "W_end_Al2", logic_w, false, 0);


	//____________________ visible attributes ____________________

    VisAtt_w = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
    VisAtt_w->SetVisibility(true);
    logic_w->SetVisAttributes(VisAtt_w);
    
    //collimator and its teeth
    VisAtt_frontcoll = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    VisAtt_frontcoll->SetVisibility(true);
    logic_frontcoll->SetVisAttributes(VisAtt_frontcoll);
    
    VisAtt_coll = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    VisAtt_coll->SetVisibility(true);
    logic_coll->SetVisAttributes(VisAtt_coll);
    
	/*
    VisAtt_coll_embed1 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //embeded heavy shielding
    VisAtt_coll_embed1->SetVisibility(true);
    logic_coll_embed1->SetVisAttributes(VisAtt_coll_embed1);
    
    VisAtt_coll_embed2 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //embeded heavy shielding
    VisAtt_coll_embed2->SetVisibility(true);
    logic_coll_embed2->SetVisAttributes(VisAtt_coll_embed2);
    
    VisAtt_coll_embed3 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //embeded heavy shielding
    VisAtt_coll_embed3->SetVisibility(true);
    logic_coll_embed3->SetVisAttributes(VisAtt_coll_embed3);
    
    VisAtt_coll_embed4= new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //embeded heavy shielding
    VisAtt_coll_embed4->SetVisibility(true);
    logic_coll_embed4->SetVisAttributes(VisAtt_coll_embed4);
    */
	VisAtt_coll_embed5 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //embeded heavy shielding
    VisAtt_coll_embed5->SetVisibility(true);
    logic_coll_embed5->SetVisAttributes(VisAtt_coll_embed5);
    /*
    VisAtt_coll_embed_Al1 = new G4VisAttributes(G4Colour(1.0,0.647,0.0)); //embeded heavy shielding
    VisAtt_coll_embed_Al1->SetVisibility(true);
    logic_coll_embed_Al1->SetVisAttributes(VisAtt_coll_embed_Al1);
    
    VisAtt_coll_embed_Al2 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));//embeded heavy shielding
    VisAtt_coll_embed_Al2->SetVisibility(true);
    logic_coll_embed_Al2->SetVisAttributes(VisAtt_coll_embed_Al2);
    
    VisAtt_coll_embed_Al3 = new G4VisAttributes(G4Colour(1.0,0.647,0.0)); //embeded heavy shielding
    VisAtt_coll_embed_Al3->SetVisibility(true);
    logic_coll_embed_Al3->SetVisAttributes(VisAtt_coll_embed_Al3);
    
    VisAtt_coll_embed_Al4 = new G4VisAttributes(G4Colour(1.0,0.647,0.0)); //embeded heavy shielding
    VisAtt_coll_embed_Al4->SetVisibility(true);
    logic_coll_embed_Al4->SetVisAttributes(VisAtt_coll_embed_Al4);
    */

    VisAtt_coll_tooth_W1 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //W1
    VisAtt_coll_tooth_W1->SetVisibility(true);
    logic_coll_tooth_W1->SetVisAttributes(VisAtt_coll_tooth_W1);
    
    VisAtt_coll_tooth_W2 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //W2
    VisAtt_coll_tooth_W2->SetVisibility(true);
    logic_coll_tooth_W2->SetVisAttributes(VisAtt_coll_tooth_W2);
    
    VisAtt_coll_tooth_W3 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //W3
    VisAtt_coll_tooth_W3->SetVisibility(true);
    logic_coll_tooth_W3->SetVisAttributes(VisAtt_coll_tooth_W3);
    
    VisAtt_coll_tooth_W4 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //W4
    VisAtt_coll_tooth_W4->SetVisibility(true);
    logic_coll_tooth_W4->SetVisAttributes(VisAtt_coll_tooth_W4);
    
    VisAtt_coll_tooth_W5 = new G4VisAttributes(G4Colour(0.0,1.0,1.0));  //W5
    VisAtt_coll_tooth_W5->SetVisibility(true);
    logic_coll_tooth_W5->SetVisAttributes(VisAtt_coll_tooth_W5);
 

    VisAtt_al_ann1 = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    VisAtt_al_ann1->SetVisibility(true);
    logic_al_ann1->SetVisAttributes(VisAtt_al_ann1);

//	VisAtt_al_ann2 = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
//    VisAtt_al_ann2->SetVisibility(true);
//    logic_al_ann2->SetVisAttributes(VisAtt_al_ann2);

    VisAtt_al_chm = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    VisAtt_al_chm->SetVisibility(true);
    logic_al_chm->SetVisAttributes(VisAtt_al_chm);

    VisAtt_al_end = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
    VisAtt_al_end->SetVisibility(true);
    logic_al_end->SetVisAttributes(VisAtt_al_end);
    
 //   VisAtt_al_endann = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
 //   VisAtt_al_endann->SetVisibility(true);
 //   logic_al_endann->SetVisAttributes(VisAtt_al_endann);

    VisAtt_be = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    VisAtt_be->SetVisibility(true);
    logic_be->SetVisAttributes(VisAtt_be);
	
    VisAtt_w_ann = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_ann->SetVisibility(true);
	logic_w_ann->SetVisAttributes(VisAtt_w_ann);
    logic_w_ann1->SetVisAttributes(VisAtt_w_ann);
	
    
/*    VisAtt_w_ann1 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_ann1->SetVisibility(true);
    logic_w_ann1->SetVisAttributes(VisAtt_w_ann1);
    
    VisAtt_w_ann2 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_ann2->SetVisibility(true);
    logic_w_ann2->SetVisAttributes(VisAtt_w_ann2);
    
    VisAtt_w_ann3 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_ann3->SetVisibility(true);
    logic_w_ann3->SetVisAttributes(VisAtt_w_ann3);
    
    VisAtt_w_ann4 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_ann4->SetVisibility(true);
    logic_w_ann4->SetVisAttributes(VisAtt_w_ann4);*/
    
    VisAtt_w_ann_Al1 = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    VisAtt_w_ann_Al1->SetVisibility(true);
    logic_w_ann_Al1->SetVisAttributes(VisAtt_w_ann_Al1);
    
    /*
	VisAtt_w_ann_Al2 = new G4VisAttributes(G4Colour(1.0, 0.647, 0.0));
    VisAtt_w_ann_Al2->SetVisibility(true);
    logic_w_ann_Al2->SetVisAttributes(VisAtt_w_ann_Al2);
    
    VisAtt_w_ann_Al3 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_ann_Al3->SetVisibility(true);
    logic_w_ann_Al3->SetVisAttributes(VisAtt_w_ann_Al3);
*/

    VisAtt_w_chm = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
    VisAtt_w_chm->SetVisibility(true);
    logic_w_chm->SetVisAttributes(VisAtt_w_chm);
    
 /*   VisAtt_w_chm1 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_chm1->SetVisibility(true);
    logic_w_chm1->SetVisAttributes(VisAtt_w_chm1);
    
    VisAtt_w_chm2 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_chm2->SetVisibility(true);
    logic_w_chm2->SetVisAttributes(VisAtt_w_chm2);
    
    VisAtt_w_chm3 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_chm3->SetVisibility(true);
    logic_w_chm3->SetVisAttributes(VisAtt_w_chm3);
    
    VisAtt_w_chm4 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_chm4->SetVisibility(true);
    logic_w_chm4->SetVisAttributes(VisAtt_w_chm4);
    
    VisAtt_w_chm5 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_chm5->SetVisibility(true);
    logic_w_chm5->SetVisAttributes(VisAtt_w_chm5);
    
    VisAtt_w_chm_Al1 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_chm_Al1->SetVisibility(true);
    logic_w_chm_Al1->SetVisAttributes(VisAtt_w_chm_Al1);
    
    VisAtt_w_chm_Al2 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_chm_Al2->SetVisibility(true);
    logic_w_chm_Al2->SetVisAttributes(VisAtt_w_chm_Al2);
    
    VisAtt_w_chm_Al3 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_chm_Al3->SetVisibility(true);
    logic_w_chm_Al3->SetVisAttributes(VisAtt_w_chm_Al3);
    
    VisAtt_w_chm_Al4 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_chm_Al4->SetVisibility(true);
    logic_w_chm_Al4->SetVisAttributes(VisAtt_w_chm_Al4);
*/
    VisAtt_w_end = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_end->SetVisibility(true);
    logic_w_end->SetVisAttributes(VisAtt_w_end);
    
 /*   VisAtt_w_end1 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_end1->SetVisibility(true);
    logic_w_end1->SetVisAttributes(VisAtt_w_end1);
    
    VisAtt_w_end2 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_end2->SetVisibility(true);
    logic_w_end2->SetVisAttributes(VisAtt_w_end2);
    
    VisAtt_w_end3 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_end3->SetVisibility(true);
    logic_w_end3->SetVisAttributes(VisAtt_w_end3);
    
    VisAtt_w_end4 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_end4->SetVisibility(true);
    logic_w_end4->SetVisAttributes(VisAtt_w_end4);
    
    VisAtt_w_end5 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_end5->SetVisibility(true);
    logic_w_end5->SetVisAttributes(VisAtt_w_end5);*/
    
    VisAtt_w_end_Al1 = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    VisAtt_w_end_Al1->SetVisibility(true);
    logic_w_end_Al1->SetVisAttributes(VisAtt_w_end_Al1);
    
    VisAtt_w_end_Al2 = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
    VisAtt_w_end_Al2->SetVisibility(true);
    logic_w_end_Al2->SetVisAttributes(VisAtt_w_end_Al2);
    
	/*
    VisAtt_w_end_Al3 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_end_Al3->SetVisibility(true);
    logic_w_end_Al3->SetVisAttributes(VisAtt_w_end_Al3);
    
    VisAtt_w_end_Al4 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_end_Al4->SetVisibility(true);
    logic_w_end_Al4->SetVisAttributes(VisAtt_w_end_Al4);
    
    //VisAtt_w_endenh = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    //VisAtt_w_endenh->SetVisibility(true);
    //logic_w_endenh->SetVisAttributes(VisAtt_w_endenh);
    
    VisAtt_w_endenh1 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_endenh1->SetVisibility(true);
    logic_w_endenh1->SetVisAttributes(VisAtt_w_endenh1);

    VisAtt_w_endenh2 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_endenh2->SetVisibility(true);
    logic_w_endenh2->SetVisAttributes(VisAtt_w_endenh2);
    
    VisAtt_w_endenh3 = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
    VisAtt_w_endenh3->SetVisibility(true);
    logic_w_endenh3->SetVisAttributes(VisAtt_w_endenh3);
    
    VisAtt_w_endenh_Al1 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_endenh_Al1->SetVisibility(true);
    logic_w_endenh_Al1->SetVisAttributes(VisAtt_w_endenh_Al1);
    
    VisAtt_w_endenh_Al2 = new G4VisAttributes(G4Colour(1.0,0.647,0.0));
    VisAtt_w_endenh_Al2->SetVisibility(true);
    logic_w_endenh_Al2->SetVisAttributes(VisAtt_w_endenh_Al2);
*/

    VisAtt_d1 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    VisAtt_d1->SetVisibility(true);
    logic_d1->SetVisAttributes(VisAtt_d1);
    

    
    VisAtt_d2 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    VisAtt_d2->SetVisibility(true);
    logic_d2->SetVisAttributes(VisAtt_d2);
    
  

    VisAtt_d3 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    VisAtt_d3->SetVisibility(true);
    logic_d3->SetVisAttributes(VisAtt_d3);
   
    
    VisAtt_d4 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
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

/*
    VisAtt_fastener1 = new G4VisAttributes(G4Colour(0.83,0.83,0.83));
    VisAtt_fastener1->SetVisibility(true);
    logic_fastener1->SetVisAttributes(VisAtt_fastener1);
    
    VisAtt_fastener2 = new G4VisAttributes(G4Colour(0.83,0.83,0.83));
    VisAtt_fastener2->SetVisibility(true);
    logic_fastener2->SetVisAttributes(VisAtt_fastener2);
    
    VisAtt_fastener3 = new G4VisAttributes(G4Colour(0.83,0.83,0.83));
    VisAtt_fastener3->SetVisibility(true);
    logic_fastener3->SetVisAttributes(VisAtt_fastener3);

	VisAtt_d5 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    VisAtt_d5->SetVisibility(true);
    logic_d5->SetVisAttributes(VisAtt_d5);
	
	VisAtt_d6 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    VisAtt_d6->SetVisibility(true);
    logic_d6->SetVisAttributes(VisAtt_d6);

	VisAtt_d7 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    VisAtt_d7->SetVisibility(true);
    logic_d7->SetVisAttributes(VisAtt_d7);

	VisAtt_d8 = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    VisAtt_d8->SetVisibility(true);
    logic_d8->SetVisAttributes(VisAtt_d8);
*/

   // designate as sensitive detectors 

  sdManager = G4SDManager::GetSDMpointer();

  //BEGIN OLD CODE FOR 2 SENSITIVE DETECTORS
 //COMMENTED BY DLT ON 3 APRIL 2008
 //G4String SiliconSensDetNames[2] = {"detector1","detector2"};
 //
 // SiSD[0] = new SetSensDet(SiliconSensDetNames[0]);
 // SiSD[1] = new SetSensDet(SiliconSensDetNames[1]);
 // sdManager->AddNewDetector(SiSD[0]);
 // sdManager->AddNewDetector(SiSD[1]);
 // logic_d1->SetSensitiveDetector(SiSD[0]);
 // logic_d2->SetSensitiveDetector(SiSD[1]);
 //END OLD CODE FOR 2 DETECTORS

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
 /*SiSD[9] = new SetSensDet(SiliconSensDetNames[9]);
 SiSD[10] = new SetSensDet(SiliconSensDetNames[10]);
 SiSD[11] = new SetSensDet(SiliconSensDetNames[11]);
 SiSD[12] = new SetSensDet(SiliconSensDetNames[12]);
 SiSD[13] = new SetSensDet(SiliconSensDetNames[13]);
 SiSD[14] = new SetSensDet(SiliconSensDetNames[14]);
 SiSD[15] = new SetSensDet(SiliconSensDetNames[15]);
 SiSD[16] = new SetSensDet(SiliconSensDetNames[16]);*/

 sdManager->AddNewDetector(SiSD[0]);
 sdManager->AddNewDetector(SiSD[1]);
 sdManager->AddNewDetector(SiSD[2]);
 sdManager->AddNewDetector(SiSD[3]);
 sdManager->AddNewDetector(SiSD[4]);
 sdManager->AddNewDetector(SiSD[5]);
 sdManager->AddNewDetector(SiSD[6]);
 sdManager->AddNewDetector(SiSD[7]);
 sdManager->AddNewDetector(SiSD[8]);
 /*sdManager->AddNewDetector(SiSD[9]);
 sdManager->AddNewDetector(SiSD[10]);
 sdManager->AddNewDetector(SiSD[11]);
 sdManager->AddNewDetector(SiSD[12]);
 sdManager->AddNewDetector(SiSD[13]);
 sdManager->AddNewDetector(SiSD[14]);
 sdManager->AddNewDetector(SiSD[15]);
 sdManager->AddNewDetector(SiSD[16]);*/

 logic_d1->SetSensitiveDetector(SiSD[0]);
 //logic_d1outer->SetSensitiveDetector(SiSD[1]);
 logic_d2->SetSensitiveDetector(SiSD[1]);
 //logic_d2outer->SetSensitiveDetector(SiSD[3]);
 logic_d3->SetSensitiveDetector(SiSD[2]);     
 //logic_d3outer->SetSensitiveDetector(SiSD[5]);
 logic_d4->SetSensitiveDetector(SiSD[3]);
 //logic_d4outer->SetSensitiveDetector(SiSD[7]);
 logic_d5->SetSensitiveDetector(SiSD[4]);
 //logic_d5outer->SetSensitiveDetector(SiSD[9]);
 logic_d6->SetSensitiveDetector(SiSD[5]);
 //logic_d6outer->SetSensitiveDetector(SiSD[11]);
 logic_d7->SetSensitiveDetector(SiSD[6]);
 //logic_d7outer->SetSensitiveDetector(SiSD[13]);
 logic_d8->SetSensitiveDetector(SiSD[7]);
 //logic_d8outer->SetSensitiveDetector(SiSD[15]);
 logic_d9->SetSensitiveDetector(SiSD[8]);
 //END NEW CODE FOR 16 DETECTORS
	return physi_w;
}

void DetectorConstruction::SetWindowMaterial(G4String materialChoice) 
{

  if (materialChoice == "Al") {
     logic_be -> SetMaterial(windowMaterialAl);
     //windowMaterial = windowMaterialAl;
     
     G4cout 
          << G4endl 
          << "----> The window is made of " << materialChoice << G4endl;
    G4cout 
          << G4endl 
          << "----> The window is made of " << logic_be->GetMaterial() << G4endl; 
    }
  if (materialChoice == "Be") {
     logic_be -> SetMaterial(windowMaterialBe);
     //windowMaterial = windowMaterialBe;
     
     G4cout 
          << G4endl 
          << "----> The window is made of " << materialChoice << G4endl; 
  }

    G4RunManager::GetRunManager()->GeometryHasBeenModified();
    //G4RunManager::GetRunManager()->ReinitializeGeometry();
    G4cout 
          << G4endl 
          << "----> The window is made of " << logic_be->GetMaterial() << G4endl; 
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
