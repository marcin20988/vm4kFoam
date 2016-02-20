/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "vm4kModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// const dataType Foam::vm4kModel::staticData();


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vm4kModel::vm4kModel(
    const twoPhaseSystem& fluid, const dimensionedVector& g):
  fluid_(fluid),
  mesh_(fluid.mesh()),
  vm4kModelDict_
  (
   IOdictionary
   (
    IOobject
    (
     "vm4kModelProperties",
     fluid_.mesh().time().constant(),
     fluid_.mesh(),
     IOobject::MUST_READ_IF_MODIFIED,
     IOobject::NO_WRITE
    )
   )
  ),  
  g_(g),
  continuousPhaseName_(vm4kModelDict_.lookup("continuousPhase")),
  dispersedPhaseName_(vm4kModelDict_.lookup("dispersedPhase")),
  T_
  (
    IOobject
    (
      "kineticTemperature",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", pow(dimLength, 2) / pow(dimTime, 2), 0.0)
  ),
  dissipation_
  (
    IOobject
    (
      "kineticDissipation",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", pow(dimLength, 2) / pow(dimTime, 3), 0.0)
  ),
  rho_
  (
    IOobject
    (
      "kineticDissipation",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", dimDensity, 0.0)
  ),
  F_
  (
    IOobject
    (
      "meanDrag",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedVector("a", dimVelocity / dimTime, vector(0.0, 0.0, 0.0))
  ),
  E1_
  (
    IOobject
    (
      "E1",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", dimless, 0.0)
  ),
  E2_
  (
    IOobject
    (
      "E2",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", pow(dimVelocity, -2), 0.0)
  ),
  cd_
  (
    IOobject
    (
      "cd",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", dimless / dimTime, 0.0)
  ),
  tau_
  (
    IOobject
    (
      "tau",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", dimTime, 0.0)
  ),
  aParameter_
  (
    IOobject
    (
      "aParameter",
      mesh_.time().timeName(),
      mesh_,
      IOobject::NO_READ,
      IOobject::AUTO_WRITE
    ),
    mesh_,
    dimensionedScalar("a", dimless, 0.0)
  )
{
  Foam::Info << "test" << Foam::endl;
}

// Foam::vm4kModel::vm4kModel(const vm4kModel&)
// {}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

// Foam::autoPtr<Foam::vm4kModel>
// Foam::vm4kModel::New()
// {
//     return autoPtr<vm4kModel>(new vm4kModel);
// }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

// Foam::vm4kModel::~vm4kModel()
// {}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::vm4kModel::update
(
    const volScalarField& T,
    const volScalarField& epsilon,
    const volScalarField& nu,
    const volVectorField& Urel,
    const volScalarField& rho
)
{
  rho_ = rho;
  T_ = T;
  dissipation_ = epsilon;
  cd_ = - fluid_.dragCoeff() 
      / dispersedPhase().rho() * dispersedPhase();
  F_ = cd_ * Urel;
  tau_ = nu / (T * (1.0f + 2.0f * cd_ * tau_));
  tau_ = max(tau_, dimensionedScalar("a", dimTime, 1e-04));
  E1_ = - dissipation_ * tau_ / T_ / (1.0f - 4.0f * cd_ * tau_);
  E2_ = - E1_ / 3.0f / T_;

  aParameter_ = dispersedPhase().d() * sqrt(F_& F_) / T_;
};

const Foam::phaseModel& Foam::vm4kModel::dispersedPhase() const
{
  if(fluid_.phase1().name() == dispersedPhaseName_)
  {
    return fluid_.phase1();
  }
  else if(fluid_.phase2().name() == dispersedPhaseName_)
  {
    return fluid_.phase2();
  }
  else
  {
    FatalErrorIn("phaseModel::dispersedPhae")
      << "there is no phase: " << dispersedPhaseName_ << endl
      << exit(FatalError);
    return fluid_.phase2();
  }
};

Foam::tmp<Foam::volScalarField> Foam::vm4kModel::kineticStress() const
{
    return rho_ * T_ * (1.0f + 7.0f * E1_ / 3.0f);
}

Foam::tmp<Foam::volScalarField> 
Foam::vm4kModel::collisionalStressScalar() const
{
    return 
      g0() * T_ 
      * 24.0f / Foam::constant::mathematical::pi
      * rho_ * dispersedPhase()
      * (
        4.0f * Foam::constant::mathematical::pi / 3.0f 
        + 2.0 * Foam::constant::mathematical::pi / 15 * pow(aParameter_, 2)
      );
}
// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

// void Foam::vm4kModel::operator=(const vm4kModel& rhs)
// {
//     // Check for assignment to self
//     if (this == &rhs)
//     {
//         FatalErrorIn("Foam::vm4kModel::operator=(const Foam::vm4kModel&)")
//             << "Attempted assignment to self"
//             << abort(FatalError);
//     }
// }

const Foam::phaseModel& Foam::vm4kModel::continuousPhase() const
{
    return fluid_.otherPhase(dispersedPhase());
}

Foam::tmp<Foam::volScalarField> Foam::vm4kModel::g0() const
{
  const volScalarField alpha = dispersedPhase();
  //return (2.0 - alpha) / (2.0 * pow(1.0 - alpha, 3));
  return 1.0 / (1.0 - alpha)
            + 3.0 * alpha / (2.0 * sqr(1.0 - alpha))
            + sqr(alpha) / (2.0 * pow3(1.0 - alpha))
            ;
}
// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //


// ************************************************************************* //
