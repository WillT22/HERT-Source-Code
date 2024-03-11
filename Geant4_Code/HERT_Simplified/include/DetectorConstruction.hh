
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
	G4LogicalVolume* logic_S1;
	G4LogicalVolume* logic_frontcoll;
	G4LogicalVolume* logic_coll;
	G4LogicalVolume* logic_coll_embed1;
	G4LogicalVolume* logic_coll_embed2;
	G4LogicalVolume* logic_coll_embed3;
	G4LogicalVolume* logic_coll_embed4;
	G4LogicalVolume* logic_coll_embed5;
	G4LogicalVolume* logic_coll_embed_Al1;
	G4LogicalVolume* logic_coll_embed_Al2;
	G4LogicalVolume* logic_coll_embed_Al3;
	G4LogicalVolume* logic_coll_embed_Al4;
	G4LogicalVolume* logic_coll_tooth_W1;
	G4LogicalVolume* logic_coll_tooth_W2;
	G4LogicalVolume* logic_coll_tooth_W3;
	G4LogicalVolume* logic_coll_tooth_W4;
	G4LogicalVolume* logic_coll_tooth_W5;
	G4LogicalVolume* logic_coll_tooth_W6;
	G4LogicalVolume* logic_coll_tooth_W7;
	G4LogicalVolume* logic_al_ann1;
	G4LogicalVolume* logic_al_chm;
	G4LogicalVolume* logic_al_end;
	G4LogicalVolume* logic_al_endann;
	G4LogicalVolume* logic_BeWin;
	G4LogicalVolume* logic_w_ann;
	G4LogicalVolume* logic_w_ann1;
	G4LogicalVolume* logic_w_ann2;
	G4LogicalVolume* logic_w_ann3;
	G4LogicalVolume* logic_w_ann4;
	G4LogicalVolume* logic_w_ann_Al1;
	G4LogicalVolume* logic_w_ann_Al2;
	G4LogicalVolume* logic_w_ann_Al3;
	G4LogicalVolume* logic_w_chm;
	G4LogicalVolume* logic_w_chm1;
	G4LogicalVolume* logic_w_chm2;
	G4LogicalVolume* logic_w_chm3;
	G4LogicalVolume* logic_w_chm4;
	G4LogicalVolume* logic_w_chm5;
	G4LogicalVolume* logic_w_chm_Al1;
	G4LogicalVolume* logic_w_chm_Al2;
	G4LogicalVolume* logic_w_chm_Al3;
	G4LogicalVolume* logic_w_chm_Al4;
	G4LogicalVolume* logic_w_end;
	G4LogicalVolume* logic_w_end1;
	G4LogicalVolume* logic_w_end2;
	G4LogicalVolume* logic_w_end3;
	G4LogicalVolume* logic_w_end4;
	G4LogicalVolume* logic_w_end5;
	G4LogicalVolume* logic_w_end_Al1;
	G4LogicalVolume* logic_w_end_Al2;
	G4LogicalVolume* logic_w_end_Al3;
	G4LogicalVolume* logic_w_end_Al4;
	G4LogicalVolume* logic_w_endenh;
	G4LogicalVolume* logic_w_endenh1;
	G4LogicalVolume* logic_w_endenh2;
	G4LogicalVolume* logic_w_endenh3;
	G4LogicalVolume* logic_w_endenh_Al1;
	G4LogicalVolume* logic_w_endenh_Al2;
	G4LogicalVolume* logic_fastener1;
	G4LogicalVolume* logic_fastener2;
	G4LogicalVolume* logic_fastener3;
	G4LogicalVolume* logic_d1;
    G4LogicalVolume* logic_d2;  
	G4LogicalVolume* logic_d3;
    G4LogicalVolume* logic_d4;
	G4LogicalVolume* logic_d5;
	G4LogicalVolume* logic_d6;
	G4LogicalVolume* logic_d7;
	G4LogicalVolume* logic_d8;
	G4LogicalVolume* logic_d9;

    // Physical volumes
	G4VPhysicalVolume* physi_w;
	G4VPhysicalVolume* physi_S1;
	G4VPhysicalVolume* physi_frontcoll;
	G4VPhysicalVolume* physi_coll;
	G4VPhysicalVolume* physi_coll_embed1;
	G4VPhysicalVolume* physi_coll_embed2;
	G4VPhysicalVolume* physi_coll_embed3;
	G4VPhysicalVolume* physi_coll_embed4;
	G4VPhysicalVolume* physi_coll_embed5;
	G4VPhysicalVolume* physi_coll_embed_Al1;
	G4VPhysicalVolume* physi_coll_embed_Al2;
	G4VPhysicalVolume* physi_coll_embed_Al3;
	G4VPhysicalVolume* physi_coll_embed_Al4;
	G4VPhysicalVolume* physi_coll_tooth_W1;
	G4VPhysicalVolume* physi_coll_tooth_W2;
	G4VPhysicalVolume* physi_coll_tooth_W3;
	G4VPhysicalVolume* physi_coll_tooth_W4;
	G4VPhysicalVolume* physi_coll_tooth_W5;
	G4VPhysicalVolume* physi_coll_tooth_W6;
	G4VPhysicalVolume* physi_coll_tooth_W7;
	G4VPhysicalVolume* physi_al_ann1;
	G4VPhysicalVolume* physi_al_chm;
	G4VPhysicalVolume* physi_al_end;
	G4VPhysicalVolume* physi_al_endann;
	G4VPhysicalVolume* physi_be;
	G4VPhysicalVolume* physi_w_ann;
	G4VPhysicalVolume* physi_w_ann1;
	G4VPhysicalVolume* physi_w_ann2;
	G4VPhysicalVolume* physi_w_ann3;
	G4VPhysicalVolume* physi_w_ann4;
	G4VPhysicalVolume* physi_w_ann_Al1;
	G4VPhysicalVolume* physi_w_ann_Al2;
	G4VPhysicalVolume* physi_w_ann_Al3;
	G4VPhysicalVolume* physi_w_chm;
	G4VPhysicalVolume* physi_w_chm1;
	G4VPhysicalVolume* physi_w_chm2;
	G4VPhysicalVolume* physi_w_chm3;
	G4VPhysicalVolume* physi_w_chm4;
	G4VPhysicalVolume* physi_w_chm5;
	G4VPhysicalVolume* physi_w_chm_Al1;
	G4VPhysicalVolume* physi_w_chm_Al2;
	G4VPhysicalVolume* physi_w_chm_Al3;
	G4VPhysicalVolume* physi_w_chm_Al4;
	G4VPhysicalVolume* physi_w_end;
	G4VPhysicalVolume* physi_w_end1;
	G4VPhysicalVolume* physi_w_end2;
	G4VPhysicalVolume* physi_w_end3;
	G4VPhysicalVolume* physi_w_end4;
	G4VPhysicalVolume* physi_w_end5;
	G4VPhysicalVolume* physi_w_end_Al1;
	G4VPhysicalVolume* physi_w_end_Al2;
	G4VPhysicalVolume* physi_w_end_Al3;
	G4VPhysicalVolume* physi_w_end_Al4;
	G4VPhysicalVolume* physi_w_endenh;
	G4VPhysicalVolume* physi_w_endenh1;
	G4VPhysicalVolume* physi_w_endenh2;
	G4VPhysicalVolume* physi_w_endenh3;
	G4VPhysicalVolume* physi_w_endenh_Al1;
	G4VPhysicalVolume* physi_w_endenh_Al2;
	G4VPhysicalVolume* physi_fastener1;
	G4VPhysicalVolume* physi_fastener2;
	G4VPhysicalVolume* physi_fastener3;

	G4VPhysicalVolume* physi_d1;
    G4VPhysicalVolume* physi_d2;
	G4VPhysicalVolume* physi_d3;
	G4VPhysicalVolume* physi_d4;
	G4VPhysicalVolume* physi_d5;
	G4VPhysicalVolume* physi_d6;
	G4VPhysicalVolume* physi_d7;
	G4VPhysicalVolume* physi_d8;
	G4VPhysicalVolume* physi_d9;
    
	// Visible attributes
	G4VisAttributes* VisAtt_w;
	G4VisAttributes* VisAtt_TaColl;
	G4VisAttributes* VisAtt_tungsten;
	G4VisAttributes* VisAtt_al;
	G4VisAttributes* VisAtt_be;
	G4VisAttributes* VisAtt_w_chm;
	G4VisAttributes* VisAtt_StStl;
	G4VisAttributes* VisAtt_d1;
    G4VisAttributes* VisAtt_d2;
	G4VisAttributes* VisAtt_d3;
	G4VisAttributes* VisAtt_d4;
	G4VisAttributes* VisAtt_d5;
	G4VisAttributes* VisAtt_d6;
	G4VisAttributes* VisAtt_d7;
	G4VisAttributes* VisAtt_d8;
	G4VisAttributes* VisAtt_d9;
	G4VisAttributes* VisAtt_coll_tooth_W1;
	G4VisAttributes* VisAtt_coll_tooth_W2;
	G4VisAttributes* VisAtt_coll_tooth_W3;
	G4VisAttributes* VisAtt_coll_tooth_W4;
	G4VisAttributes* VisAtt_coll_tooth_W5;
	G4VisAttributes* VisAtt_coll_embed5;
	G4VisAttributes* VisAtt_coll;
	G4VisAttributes* VisAtt_frontcoll;
	G4VisAttributes* VisAtt_al_ann1;
	G4VisAttributes* VisAtt_al_chm;
	G4VisAttributes* VisAtt_al_end;

	G4VisAttributes* VisAtt_w_ann;
	G4VisAttributes* VisAtt_w_end;
	G4VisAttributes* VisAtt_w_ann_Al1;
	G4VisAttributes* VisAtt_w_end_Al1;
	G4VisAttributes* VisAtt_w_end_Al2;

	
	 G4SDManager* sdManager;
public:	 
     G4VSensitiveDetector* SiSD[9]; 
    
};

#endif

