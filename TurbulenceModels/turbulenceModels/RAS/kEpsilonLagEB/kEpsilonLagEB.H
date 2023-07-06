/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::RASModels::kEpsilonLagEB

Group
    grpRASTurbulence

Description
    The k-epsilon Lag Elliptic-Blending turbulence closure model for incompressible and
    compressible flows.

    The model is a three-transport-equation linear-eddy-viscosity turbulence
    closure model alongside an elliptic blending equation.

    \heading Input fields
    \plaintable
        k         | Turbulent kinetic energy [m2/s2]
        epsilon   | Turbulent kinetic energy dissipation rate [m2/s3]
        phit      | Lag between stress-strain [-]
        f         | Elliptic blending factor [-]
    \endplaintable

    Reference:
    \verbatim
        Standard model (Tag:TLLP):
            Laurence
            Tunstall, R., Lardeau, S., Laurence, D., Prosser, R. (2016).
            An Elliptic Blending Lag Model for Flows in Thermal-Hydraulics Systems
            ETMM11 Conference
            https://www.researchgate.net/publication/308413285_AN_ELLIPTIC_BLENDING_LAG_MODEL_FOR_FLOWS_IN_THERMAL-HYDRAULICS_SYSTEMS
    \endverbatim

    The default model coefficients are:
    \verbatim
        kEpsilonLagEBCoeffs
        {
            Cmu         0.22;       // Turbulent viscosity constant
            C1          1.7;
            C3          0.8;
            C4          0.625;
            C5          0.2;
            C1s         0.9;
            C2s         0.65;       
            CK          2.3;        // Model constant for E
            Ct          4.0;        // Model constant for T
            Ceps1       1.44;       // Model constant for epsilon
            Ceps2       1.9;        // Model constant for epsilon
            CL          0.164;      // Model constant for L
            Ceta        75.0;       // Model constant for L
            sigmaK      1.0;        // Turbulent Prandtl number for k
            sigmaEps    1.2;        // Turbulent Prandtl number for epsilon
            sigmaPhit   1.0;        // Turbulent Prandtl number for phit
        }
    \endverbatim

Note
    The name of the variable representing the lag between stress and strain is 'phi'.
    However, the name 'phi' preexisted in OpenFOAM; therefore, this name was
    replaced by 'phit' herein.

SourceFiles
    kEpsilonLagEB.C

SeeAlso
    kEpsilon.C

\*---------------------------------------------------------------------------*/

#ifndef kEpsilonLagEB_H
#define kEpsilonLagEB_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                        Class kEpsilonLagEB Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kEpsilonLagEB
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        kEpsilonLagEB(const kEpsilonLagEB&) = delete;

        //- No copy assignment
        void operator=(const kEpsilonLagEB&) = delete;


protected:

    // Protected Data
        // Model coefficients

            dimensionedScalar Cmu_;
            dimensionedScalar C1_;
            dimensionedScalar C3_;
            dimensionedScalar C4_;
            dimensionedScalar C5_;
            dimensionedScalar C1s_;
            dimensionedScalar C3s_;
            dimensionedScalar CK_;
            dimensionedScalar Ct_;
            dimensionedScalar Ceps1_;
            dimensionedScalar Ceps2_;
            dimensionedScalar CL_;
            dimensionedScalar Ceta_;
            dimensionedScalar sigmaK_;
            dimensionedScalar sigmaEps_;
            dimensionedScalar sigmaPhit_;


        // Fields

            //- Turbulent kinetic energy [m2/s2]
            volScalarField k_;

            //- Turbulent kinetic energy dissipation rate [m2/s3]
            volScalarField epsilon_;

            //- Lag between stress-strain [-]
            volScalarField phit_;

            //- Elliptic blending factor [-]
            volScalarField ebf_;

            //- Turbulent time scale [s]
            volScalarField T_;


        // Bounding values

            dimensionedScalar phitMin_;
            dimensionedScalar ebfMin_;
            dimensionedScalar TMin_;
            dimensionedScalar L2Min_;

    // Protected Member Functions

        //- Update nut with the latest available k, phit, and T
        virtual void correctNut();

        //- Return the turbulent time scale, T
        tmp<volScalarField> Ts() const;

        //- Return the turbulent length scale, L
        tmp<volScalarField> Ls() const;


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kEpsilonLagEB");


    // Constructors

        //- Construct from components
        kEpsilonLagEB
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kEpsilonLagEB() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for k (TLLP:Eq.5)
        tmp<volScalarField> DkEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DkEff",
                    this->nut_/sigmaK_ + this->nu()/2
                )
            );
        }

        //- Return the effective diffusivity for epsilon (TLLP:Eq.6)
        tmp<volScalarField> DepsilonEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DepsilonEff",
                    this->nut_/sigmaEps_ + this->nu()/2
                )
            );
        }

        //- Return the effective diffusivity for phit (TLLP:Eq.7)
        tmp<volScalarField> DphitEff() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DphitEff",
                    this->nut_/sigmaPhit_ + this->nu()/2
                )
            );
        }

        //- Return the turbulent kinetic energy field
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulent kinetic energy dissipation rate field
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the normalised wall-normal fluctuating velocity scale field
        virtual tmp<volScalarField> phit() const
        {
            return phit_;
        }

        //- Return the elliptic relaxation factor field
        virtual tmp<volScalarField> ebf() const
        {
            return ebf_;
        }

        //- Solve the transport equations and correct the turbulent viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kEpsilonLagEB.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
