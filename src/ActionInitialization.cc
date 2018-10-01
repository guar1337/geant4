#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

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
	SetUserAction(new PrimaryGeneratorAction);
}