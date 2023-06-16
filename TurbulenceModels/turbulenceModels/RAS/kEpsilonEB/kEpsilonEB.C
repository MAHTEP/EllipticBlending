/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "kEpsilonEB.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void kEpsilonEB<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*phit_*k_*T_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonEB<BasicTurbulenceModel>::Ts() const
{
    return
        sqrt
        (   
            sqr(k_/epsilon_) + sqr(Ct_)*max(this->nu()/epsilon_,
                dimensionedScalar(sqr(dimTime), Zero) )

        );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonEB<BasicTurbulenceModel>::Ls() const
{
    return
       CL_*sqrt(max
       (
        pow3(k_)/sqr(epsilon_) + sqr(Ceta_)*sqrt(max
        (
            pow3(this->nu())/epsilon_,
            dimensionedScalar(pow(dimLength,4), Zero)
        )),
        dimensionedScalar(sqr(dimLength), Zero)
       ));
       
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonEB<BasicTurbulenceModel>::kEpsilonEB
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    includeNu_
    (
        Switch::getOrAddToDict
        (
            "includeNu",
            this->coeffDict_,
            true
        )
    ),
    Cmu_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.22
        )
    ),
    C1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C1_",
            this->coeffDict_,
            1.7
        )
    ),
    C2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            0.9
        )
    ),
    CK_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CK",
            this->coeffDict_,
            2.3
        )
    ),
    Ct_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ct",
            this->coeffDict_,
            4.0
        )
    ),
    Ceps1_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps1",
            this->coeffDict_,
            1.44
        )
    ),
    Ceps2_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceps2",
            this->coeffDict_,
            1.83
        )
    ),
    CL_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CL",
            this->coeffDict_,
            0.164
        )
    ),
    Ceta_
    (
        dimensionedScalar::getOrAddToDict
        (
            "Ceta",
            this->coeffDict_,
            75.0
        )
    ),
    CT_
    (
        dimensionedScalar::getOrAddToDict
        (
            "CT",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaK_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaK",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.5
        )
    ),
    sigmaPhit_
    (
        dimensionedScalar::getOrAddToDict
        (
            "sigmaPhit",
            this->coeffDict_,
            1.0
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    phit_
    (
        IOobject
        (
            IOobject::groupName("phit", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    ebf_
    (
        IOobject
        (
            IOobject::groupName("ebf", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    T_
    (
        IOobject
        (
            "T",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        this->mesh_,
        dimensionedScalar(dimTime, Zero)
    ),

    phitMin_(dimensionedScalar("phitMin", phit_.dimensions(), 1e-10)),
    ebfMin_(dimensionedScalar("ebfMin", ebf_.dimensions(), Zero)),
    TMin_(dimensionedScalar("TMin", dimTime, SMALL)),
    L2Min_(dimensionedScalar("L2Min", sqr(dimLength), SMALL))
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);
    bound(phit_, phitMin_);
    bound(ebf_, ebfMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }

    if
    (
        mag(sigmaK_.value()) < VSMALL
     || mag(sigmaEps_.value()) < VSMALL
     || mag(sigmaPhit_.value()) < VSMALL
    )
    {
        FatalErrorInFunction
            << "Non-zero values are required for the model constants:" << nl
            << "sigmaK = " << sigmaK_ << nl
            << "sigmaEps = " << sigmaEps_ << nl
            << "sigmaPhit = " << sigmaPhit_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kEpsilonEB<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        includeNu_.readIfPresent("includeNu", this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        CK_.readIfPresent(this->coeffDict());
        Ct_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        CT_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        sigmaPhit_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kEpsilonEB<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Construct local convenience references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))
    );

    tmp<volSymmTensorField> tS(symm(fvc::grad(U)));
    volScalarField G(this->GName(), nut*(2.0*(dev(tS()) && tS())));
    tS.clear();

    T_ = Ts();
    bound(T_, TMin_);
    // compute strain, vorticity and anisotropy tensors
    tmp<volTensorField> tgradU = fvc::grad(U);
    volTensorField S
    (
        0.5*(tgradU() + tgradU().T())
    );
    tmp<volVectorField> gradf(fvc::grad(ebf_));
    volVectorField n 
    (
        fvc::grad(ebf_)/max(mag(fvc::grad(ebf_)), dimensionedScalar(dimless/dimLength, SMALL))
    );
    gradf.clear();
    volScalarField E
    (
        CK_*pow3(1-ebf_)*this->nu()*nut*sqr(fvc::div(mag(2*S&n)*n))
    );

    const volScalarField L2(type() + "L2", sqr(Ls()) + L2Min_);


    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Turbulent kinetic energy dissipation rate equation 
    // k/T ~ epsilon
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
      ==
        alpha()*rho()*Ceps1_*G()/k_()*epsilon_()
      - fvm::SuSp
        (
            (2.0/3.0*Ceps1_)*(alpha()*rho()*divU),
            epsilon_
        )
      - fvm::Sp(alpha()*rho()*Ceps2_/k_()*epsilon_(), epsilon_)
      + alpha()*rho()*E()
      + fvOptions(alpha, rho, epsilon_)
    );
    
    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation
    // epsilon/k ~ 1/Ts
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
      ==
        alpha()*rho()*G()
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*(1.0/T_()), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


// Elliptic blending function equation
    tmp<fvScalarMatrix> ebfEqn
   (
    -fvm::laplacian(ebf_)
    ==
    -fvm::SuSp(1.0/L2(), ebf_) + 1.0/L2()
   );

    ebfEqn.ref().relax();
    fvOptions.constrain(ebfEqn.ref());
    solve(ebfEqn);
    fvOptions.correct(ebf_);
    bound(ebf_, this->ebfMin_);

// Reduced stress function
    tmp<fvScalarMatrix> phitEqn
    (
        fvm::ddt(alpha, rho, phit_)
      + fvm::div(alphaRhoPhi, phit_)
      - fvm::laplacian(alpha*rho*DphitEff(), phit_)
      ==
      - fvm::SuSp
        (
            alpha()*rho()*
            (
                G()/k_()
              - (2.0/3.0)*divU
              + 0.5*(1-pow3(ebf_()))/k_()*epsilon_()
              + pow3(ebf_())/k_()*epsilon_()*(C1_ - 1.0 + C2_*(G()-(2.0/3.0)*divU*k_())/epsilon_())
            )
          , phit_
        )
        + alpha()*rho()*2.0/3.0*pow3(ebf_())/k_()*epsilon_()*(C1_ - 1.0 + C2_*G()/epsilon_())
      + fvOptions(alpha, rho, phit_)
    );

    phitEqn.ref().relax();
    fvOptions.constrain(phitEqn.ref());
    solve(phitEqn);
    fvOptions.correct(phit_);
    bound(phit_, phitMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
