/*
 * The MIT License
 *
 * Copyright (c) 1997-2019 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */


// SingleVel.cc
// One of the derived Contact classes.  This particular
// class contains methods for recapturing single velocity
// field behavior from objects belonging to multiple velocity
// fields.  The main purpose of this type of contact is to
// ensure that one can get the same answer using prescribed
// contact as can be gotten using "automatic" contact.
#include <CCA/Components/MPM/Materials/Contact/SingleVelLiquidContact.h>
#include <Core/Geometry/Vector.h>
#include <Core/Geometry/IntVector.h>
#include <Core/Grid/Variables/NCVariable.h>
#include <Core/Grid/Patch.h>
#include <Core/Grid/Level.h>
#include <Core/Grid/Variables/NodeIterator.h>
#include <Core/Grid/MaterialManager.h>
#include <Core/Grid/MaterialManagerP.h>
#include <Core/Grid/Task.h>
#include <Core/Grid/Variables/VarTypes.h>
#include <CCA/Components/MPM/Core/MPMLabel.h>
#include <CCA/Components/MPM/Core/DOUBLEMPMLabel.h>
#include <CCA/Ports/DataWarehouse.h>
#include <vector>

using namespace std;
using namespace Uintah;
using std::vector;

SingleVelLiquidContact::SingleVelLiquidContact(const ProcessorGroup* myworld,
                                   ProblemSpecP& ps, MaterialManagerP& d_sS, 
                                   MPMLabel* Mlb,MPMFlags* MFlag)
  : Contact(myworld, Mlb, MFlag, ps)
{
  // Constructor
  d_materialManager = d_sS;
  lb = Mlb;
  flag = MFlag;
  d_oneOrTwoStep = 2;
  ps->get("OneOrTwoStep",     d_oneOrTwoStep);
}

SingleVelLiquidContact::~SingleVelLiquidContact()
{
}

void SingleVelLiquidContact::outputProblemSpec(ProblemSpecP& ps)
{
  ProblemSpecP contact_ps = ps->appendChild("contact");
  contact_ps->appendElement("type","single_velocity");
  d_matls.outputProblemSpec(contact_ps);
}

void SingleVelLiquidContact::exMomInterpolated(const ProcessorGroup*,
                                         const PatchSubset* patches,
                                         const MaterialSubset* matls,
                                         DataWarehouse*,
                                         DataWarehouse* new_dw)
{
  if(d_oneOrTwoStep==2){
   string interp_type = flag->d_interpolator_type;
   int numMatls = d_materialManager->getNumMatls( "MPM" );
   ASSERTEQ(numMatls, matls->size());
   for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    Vector centerOfMassVelocity(0.0,0.0,0.0);

    // Retrieve necessary data from DataWarehouse
    std::vector<constNCVariable<double> > gmass(numMatls);
    //std::vector<NCVariable<Vector> > gvelocity(numMatls);
	std::vector<NCVariable<Vector> > gvelocityLiquid(numMatls);
    for(int m=0;m<matls->size();m++){
      int dwindex = matls->get(m);
      new_dw->get(gmass[m], lb->gMassLabel,    dwindex, patch,Ghost::None,0);
      //new_dw->getModifiable(gvelocity[m], lb->gVelocityLabel,dwindex, patch);
	  new_dw->getModifiable(gvelocityLiquid[m], lb->gVelocityLiquidLabel, dwindex, patch);
    }

    for(NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++){
      IntVector c = *iter;

      Vector centerOfMassMom(0,0,0);
      double centerOfMassMass=0.0;

      for(int n = 0; n < numMatls; n++){
        if(d_matls.requested(n)) {
          centerOfMassMom+= gvelocityLiquid[n][c] * gmass[n][c];
          centerOfMassMass+=gmass[n][c]; 
        }
      }

      // Set each field's velocity equal to the center of mass velocity
      centerOfMassVelocity=centerOfMassMom/centerOfMassMass;
      for(int n = 0; n < numMatls; n++) {
        if(d_matls.requested(n)) {
			gvelocityLiquid[n][c] = centerOfMassVelocity;
        }
      }
    }
  }
 }
}

void SingleVelLiquidContact::exMomIntegrated(const ProcessorGroup*,
                                       const PatchSubset* patches,
                                       const MaterialSubset* matls,
                                       DataWarehouse* old_dw,
                                       DataWarehouse* new_dw)
{
  int numMatls = d_materialManager->getNumMatls( "MPM" );
  ASSERTEQ(numMatls, matls->size());

  for(int p=0;p<patches->size();p++){
    const Patch* patch = patches->get(p);
    Vector zero(0.0,0.0,0.0);
    Vector centerOfMassVelocity(0.0,0.0,0.0);
    Vector centerOfMassMom(0.0,0.0,0.0);
    Vector Dvdt;
    double centerOfMassMass;

    // Retrieve necessary data from DataWarehouse
    std::vector<constNCVariable<double> > gmass(numMatls);
    std::vector<NCVariable<Vector> > gvelocity_star(numMatls);
	std::vector<NCVariable<Vector> > gvelocityLiquid_star(numMatls);
    for(int m=0;m<matls->size();m++){
     int dwi = matls->get(m);
     new_dw->get(gmass[m],lb->gMassLabel, dwi, patch, Ghost::None, 0);
     //new_dw->getModifiable(gvelocity_star[m],lb->gVelocityStarLabel, dwi,patch);
	 new_dw->getModifiable(gvelocityLiquid_star[m], double_lb->gVelocityStarLiquidLabel,
		 dwi, patch);
    }

    delt_vartype delT;
    old_dw->get(delT, lb->delTLabel, getLevel(patches));
    
    for(NodeIterator iter = patch->getNodeIterator(); !iter.done(); iter++){
      IntVector c = *iter;

      centerOfMassMom=zero;
      centerOfMassMass=0.0; 
      for(int  n = 0; n < numMatls; n++){
        if(d_matls.requested(n)) {
          centerOfMassMom+= gvelocityLiquid_star[n][c] * gmass[n][c];
          centerOfMassMass+=gmass[n][c]; 
        }
      }

      // Set each field's velocity equal to the center of mass velocity
      centerOfMassVelocity=centerOfMassMom/centerOfMassMass;
      for(int  n = 0; n < numMatls; n++){
        if(d_matls.requested(n)) {
          Dvdt = (centerOfMassVelocity - gvelocityLiquid_star[n][c])/delT;
		  gvelocityLiquid_star[n][c] = centerOfMassVelocity;
        }
      }
    }
  }
}

void SingleVelLiquidContact::addComputesAndRequiresInterpolated(SchedulerP & sched,
                                                  const PatchSet* patches,
                                                  const MaterialSet* ms)
{
  Task * t = scinew Task("SingleVelLiquidContact::exMomInterpolated", 
                      this, &SingleVelLiquidContact::exMomInterpolated);
  
  const MaterialSubset* mss = ms->getUnion();
  t->requires( Task::NewDW, lb->gMassLabel,          Ghost::None);

  //t->modifies(lb->gVelocityLabel, mss);
  t->modifies(lb->gVelocityLiquidLabel, mss);

  sched->addTask(t, patches, ms);
}

void SingleVelLiquidContact::addComputesAndRequiresIntegrated(SchedulerP & sched,
                                             const PatchSet* patches,
                                             const MaterialSet* ms) 
{
  Task * t = scinew Task("SingleVelLiquidContact::exMomIntegrated", 
                      this, &SingleVelLiquidContact::exMomIntegrated);
  
  const MaterialSubset* mss = ms->getUnion();
  t->requires(Task::OldDW, lb->delTLabel);    
  t->requires(Task::NewDW, lb->gMassLabel,              Ghost::None);

  //t->modifies(lb->gVelocityStarLabel, mss);
  t->modifies(double_lb->gVelocityStarLiquidLabel, mss);
  sched->addTask(t, patches, ms);
}
