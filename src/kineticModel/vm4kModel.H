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

Class
    Foam::vm4kModel

Description

SourceFiles
    vm4kModelI.H
    vm4kModel.C
    vm4kModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef vm4kModel_H
#define vm4kModel_H

#include "twoPhaseSystem.H"
#include "phaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


// Forward declaration of friend functions and operators
class vm4kModel;

/*---------------------------------------------------------------------------*\
                         Class vm4kModel Declaration
\*---------------------------------------------------------------------------*/

class vm4kModel
{
    // Private data
    private:
        const twoPhaseSystem& fluid_;
        const fvMesh& mesh_;
    protected: 
        const dictionary vm4kModelDict_;
    private:
        const dimensionedVector g_;
        const word continuousPhaseName_;
        const word dispersedPhaseName_;

        volScalarField T_;
        volScalarField dissipation_;
        volScalarField rho_;
        volVectorField F_;
        volScalarField E1_;
        volScalarField E2_;
        volScalarField cd_;
        volScalarField tau_;
        volScalarField aParameter_;

        //- Description of data_


    // Private Member Functions

        //- Disallow default bitwise copy construct
        // vm4kModel(const vm4kModel&);

        //- Disallow default bitwise assignment
        // void operator=(const vm4kModel&);

public:

    // Static data members

        //- Static data staticData

    // Constructors

        //- Construct null
        // vm4kModel();

        vm4kModel(const twoPhaseSystem& fluid, const dimensionedVector& g);
        

        void update(
            const volScalarField& T,
            const volScalarField& epsilon,
            const volScalarField& nu,
            const volVectorField& F,
            const volScalarField& rho
        );

        const phaseModel& dispersedPhase() const;
        const phaseModel& continuousPhase() const;

        tmp<volScalarField> kineticStress() const;
        tmp<volScalarField> g0() const;
        tmp<volScalarField> collisionalStressScalar() const;

        //- Construct from components
        // vm4kModel(const dataType& data);

        //- Construct from Istream

        //- Construct as copy
        // vm4kModel(const vm4kModel&);


    // Selectors

        //- Select null constructed
        // static autoPtr<vm4kModel> New();


    //- Destructor
    // ~vm4kModel();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        // void operator=(const vm4kModel&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
