## Information

* The clustered_dna_damage application was developed using a modified version of Topas 3.3.p1.
* Modifications were required to enable efficient geometry construction of nuclear DNA.
* The required methods were introduced in Topas 3.6.
* As of this time, Topas-nBio 1.0-beta.1 is incompatible with Topas 3.6.
* Thus, please add the following code for 2 versions of the `CreatePhysicalVolume()` to `TsVGeometryComponent.hh`.
* `TsVGeometryComponent.hh` can be found in the `include/` folder of your local Topas installation.
* Once the code has been inserted, recompile your Topas installation as per the Topas instructions.

## Code

### Include statements

* Add this code alongside the other `include` statements (near the top of the file).

```c++
//--------------------------------------------------------------------------------------------------
// Additions required for clustered_dna_damage application
#include "G4UIcommand.hh"
#include "G4PVPlacement.hh"
//--------------------------------------------------------------------------------------------------
```

### Additional methods

* Add this code into the `protected` methods section.
* E.g. between the existing `CreatePhysicalVolume` methods and the `RegisterRotation` method.

```c++
//----------------------------------------------------------------------------------------------
// Additions required for clustered_dna_damage application
G4VPhysicalVolume* CreatePhysicalVolume(const char* subComponentName, G4int copy, G4bool reuseLogical, G4LogicalVolume* lVol,
                                                              G4RotationMatrix* rot, G4ThreeVector* trans, G4LogicalVolume* parent) {
    G4String nameString = subComponentName;
    return CreatePhysicalVolume(nameString, copy, reuseLogical, lVol, rot, trans, parent);
}


G4VPhysicalVolume* CreatePhysicalVolume(G4String& subComponentName, G4int copy, G4bool reuseLogical, G4LogicalVolume* lVol,
                                                              G4RotationMatrix* rot, G4ThreeVector* trans, G4LogicalVolume* parent) {
    G4String pVolName = fName + fCopyId;
    if (subComponentName!="") pVolName += "/" + subComponentName;
    if (copy!=-1) pVolName += G4UIcommand::ConvertToString(copy);

    G4int replica;
    if (reuseLogical) replica = copy;
    else replica = 0;

    G4VPhysicalVolume* pvol = new G4PVPlacement(rot, *trans, lVol, pVolName, parent, false, replica, false);

    CheckForOverlaps(pvol);

    fRotations.push_back(rot);
    fTranslations.push_back(trans);
    fPhysicalVolumes.push_back(pvol);
    return pvol;
}
//----------------------------------------------------------------------------------------------
```