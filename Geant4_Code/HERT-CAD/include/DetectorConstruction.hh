
#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1


class G4Material;
class G4UserLimits;
class G4GlobalMagFieldMessenger;
class DetectorMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VisAttributes;

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4RunManager.hh"
#include "DetectorMessenger.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
	
    DetectorConstruction();
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();
    // Set methods
    void SetWindowMaterial (G4String );
    void SetWindowDepth  (G4double val);
  private:
        DetectorMessenger*  fDetectorMessenger;
    // Logical volumes
    G4double         windowDepth;
    G4Material*      windowMaterialAl;
    G4Material*      windowMaterialBe;
    G4Material*      windowMaterial;
	G4Material*		 AlAlloy;
    
	G4LogicalVolume* logic_w;
	G4LogicalVolume* logic_alHousing;
	G4LogicalVolume* logic_BeWin;
	G4LogicalVolume* logic_TaTooth5;
	G4LogicalVolume* logic_TaTooth4;
	G4LogicalVolume* logic_TaSpac4;
	G4LogicalVolume* logic_TaTooth3;
	G4LogicalVolume* logic_TaSpac3;
	G4LogicalVolume* logic_TaTooth2;
	G4LogicalVolume* logic_TaSpac2;
	G4LogicalVolume* logic_TaTooth1;
	G4LogicalVolume* logic_TaSpac1;
	G4LogicalVolume* logic_WFr1;
	G4LogicalVolume* logic_AlFrShim;
	G4LogicalVolume* logic_WFrIn1;
	G4LogicalVolume* logic_AlChShim1;
	G4LogicalVolume* logic_BackW;
	G4LogicalVolume* logic_EpoxCham;
	G4LogicalVolume* logic_AlBkShim1;
	G4LogicalVolume* logic_AlBkPlate;
	G4LogicalVolume* logic_AlignPin1;
	G4LogicalVolume* logic_AlignPin2;
	G4LogicalVolume* logic_AlignPin3; 
	
	G4LogicalVolume* logic_d1;
	G4LogicalVolume* logic_d2;
	G4LogicalVolume* logic_d3;
	G4LogicalVolume* logic_d4;
	G4LogicalVolume* logic_d5;
	G4LogicalVolume* logic_d6;
	G4LogicalVolume* logic_d7;
	G4LogicalVolume* logic_d8;
	G4LogicalVolume* logic_d9;
    
	G4LogicalVolume* logic_S1;

    // Physical volumes
	G4VPhysicalVolume* physi_w;

	G4VPhysicalVolume* physi_alHousing;
	G4VPhysicalVolume* physi_BeWin;
	G4VPhysicalVolume* physi_TaTooth5;
	G4VPhysicalVolume* physi_TaTooth4;
	G4VPhysicalVolume* physi_TaSpac4;
	G4VPhysicalVolume* physi_TaTooth3;
	G4VPhysicalVolume* physi_TaSpac3;
	G4VPhysicalVolume* physi_TaTooth2;
	G4VPhysicalVolume* physi_TaSpac2;
	G4VPhysicalVolume* physi_TaTooth1;
	G4VPhysicalVolume* physi_TaSpac1;
	G4VPhysicalVolume* physi_WFr1;
	G4VPhysicalVolume* physi_AlFrShim;
	G4VPhysicalVolume* physi_AlChShim1;
	G4VPhysicalVolume* physi_BackW;
	G4VPhysicalVolume* physi_EpoxCham;
	G4VPhysicalVolume* physi_AlBkShim1;
	G4VPhysicalVolume* physi_AlBkPlate;
	G4VPhysicalVolume* physi_AlignPin1;
	G4VPhysicalVolume* physi_AlignPin2;
	G4VPhysicalVolume* physi_AlignPin3;
	G4VPhysicalVolume* physi_WFrIn1;
	
	G4VPhysicalVolume* physi_d1;
	G4VPhysicalVolume* physi_d2;
	G4VPhysicalVolume* physi_d3;
	G4VPhysicalVolume* physi_d4;
	G4VPhysicalVolume* physi_d5;
	G4VPhysicalVolume* physi_d6;
	G4VPhysicalVolume* physi_d7;
	G4VPhysicalVolume* physi_d8;
	G4VPhysicalVolume* physi_d9;

	G4VPhysicalVolume* physi_S1;
	
	
    
	// Visible attributes
	G4VisAttributes* VisAtt_w;
	G4VisAttributes* VisAtt_detectors;
	G4VisAttributes* VisAtt_TaColl;
	G4VisAttributes* VisAtt_tungsten;
	G4VisAttributes* VisAtt_al;
	G4VisAttributes* VisAtt_be;
	G4VisAttributes* VisAtt_w_chm;
	G4VisAttributes* VisAtt_StStl;
	G4VisAttributes* VisAtt_fastener1;
	G4VisAttributes* VisAtt_fastener2;
	G4VisAttributes* VisAtt_fastener3;
	
	 G4SDManager* sdManager;
public:	 
     G4VSensitiveDetector* SiSD[9]; 
    
};

#endif

