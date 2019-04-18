#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
<<<<<<< HEAD
=======
#include "EventAction.hh"
#include "DetectorConstruction.hh"
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217

ActionInitialization::ActionInitialization()
 : G4VUserActionInitialization()
{
}

ActionInitialization::~ActionInitialization()
{}

void ActionInitialization::BuildForMaster() const
{
}

void ActionInitialization::Build() const
{
<<<<<<< HEAD
	SetUserAction(new PrimaryGeneratorAction{});
=======
	SetUserAction(new PrimaryGeneratorAction);
>>>>>>> a51287d9dc94593d291b50567b67b4de15ec5217
}