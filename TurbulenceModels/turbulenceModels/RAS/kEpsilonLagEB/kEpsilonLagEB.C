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

#include "kEpsilonLagEB.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

// Compute the turbulent viscosity (TLLP:Eq.10)
template<class BasicTurbulenceModel>
void kEpsilonLagEB<BasicTurbulenceModel>::correctNut()
{   
    volScalarField magS =
        sqrt(2.0)*mag(symm(fvc::grad(this->U_)));

    this->nut_ = Cmu_*phit_*k_*
        min
        (
            T_, 
            1.0
            /
            max
            (
                dimensionedScalar(pow(dimTime,-1),VSMALL),
                Cmu_.value()*sqrt(3.0)*phit_*magS
                    
            )
        );
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}

// Compute the limit turbulent time scale (TLLP:Eq.14)
template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonLagEB<BasicTurbulenceModel>::Ts() const
{
    return
        sqrt
        ( 
            sqr(k_/epsilon_) + sqr(Ct_)*max(this->nu()/epsilon_,
                dimensionedScalar(sqr(dimTime), Zero) )
        );
}

// Compute the turbulent length scale (TLLP:Eq.12)
template<class BasicTurbulenceModel>
tmp<volScalarField> kEpsilonLagEB<BasicTurbulenceModel>::Ls() const
{
    return
       CL_*sqrt(
        max(pow3(k_)/sqr(epsilon_), dimensionedScalar(sqr(dimLength), Zero))
        + sqr(Ceta_)
        *sqrt(
            max(
                pow3(this->nu())/epsilon_,
                dimensionedScalar(pow(dimLength,4), Zero)
                )
            )
       );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kEpsilonLagEB<BasicTurbulenceModel>::kEpsilonLagEB
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
            "C1",
            this->coeffDict_,
            1.7
        )
    ),
    C3_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0.8
        )
    ),
    C4_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C4",
            this->coeffDict_,
            0.625
        )
    ),
    C5_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C5",
            this->coeffDict_,
            0.2
        )
    ),
    C1s_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C1s",
            this->coeffDict_,
            0.9
        )
    ),
    C3s_
    (
        dimensionedScalar::getOrAddToDict
        (
            "C3s",
            this->coeffDict_,
            0.65
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
            "Ceps1_",
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
            1.9
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
            1.2
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

    phitMin_(dimensionedScalar("phitMin", phit_.dimensions(), VSMALL)),
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
bool kEpsilonLagEB<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        C4_.readIfPresent(this->coeffDict());
        C5_.readIfPresent(this->coeffDict());
        C1s_.readIfPresent(this->coeffDict());
        C3s_.readIfPresent(this->coeffDict());
        CK_.readIfPresent(this->coeffDict());
        Ct_.readIfPresent(this->coeffDict());
        Ceps1_.readIfPresent(this->coeffDict());
        Ceps2_.readIfPresent(this->coeffDict());
        CL_.readIfPresent(this->coeffDict());
        Ceta_.readIfPresent(this->coeffDict());
        sigmaK_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        sigmaPhit_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void kEpsilonLagEB<BasicTurbulenceModel>::correct()
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

    const volScalarField L2(type() + "L2", sqr(Ls()) + L2Min_);
 
    // Compute strain, vorticity and anisotropy tensors
    tmp<volTensorField> tgradU = fvc::grad(U);
    
    // Mean strain rate tensor
    const volTensorField S
    (
        0.5*(tgradU() + tgradU().T())
    );

    // Mean vorticity tensor
    const volTensorField W
    (
        0.5*(tgradU() - tgradU().T())
    );
    tgradU.clear();

    // Definition of components of the mean strain rate
    volScalarField Sxx_(S.component(tensor::XX));
    volScalarField Sxy_(S.component(tensor::XY));
    volScalarField Sxz_(S.component(tensor::XZ));
    volScalarField Syx_(S.component(tensor::YX));
    volScalarField Syy_(S.component(tensor::YY));
    volScalarField Syz_(S.component(tensor::YZ));
    volScalarField Szx_(S.component(tensor::ZX));
    volScalarField Szy_(S.component(tensor::ZY));
    volScalarField Szz_(S.component(tensor::ZZ));
    
    // Defition of 1/S^2*(S*DS/Dt)
   volTensorField SDS(S);

    // SDS.component(tensor::XY) = (Sxx_*(fvc::ddt(Syx_)+fvc::div(this->phi(),Syx_))
    //                           + Sxy_*(fvc::ddt(Syy_)+fvc::div(this->phi(),Syy_))
    //                           + Sxz_*(fvc::ddt(Syz_)+fvc::div(this->phi(),Syz_)))/(2.0*magSqr(S));
    // SDS.component(tensor::XZ) = (Sxx_*(fvc::ddt(Szx_)+fvc::div(this->phi(),Szx_))
    //                           + Sxy_*(fvc::ddt(Szy_)+fvc::div(this->phi(),Szy_))
    //                           + Sxz_*(fvc::ddt(Szz_)+fvc::div(this->phi(),Szz_)))/(2.0*magSqr(S));                          
    // SDS.component(tensor::YX) = (Syx_*(fvc::ddt(Sxx_)+fvc::div(this->phi(),Sxx_))
    //                           + Syy_*(fvc::ddt(Sxy_)+fvc::div(this->phi(),Sxy_))
    //                           + Syz_*(fvc::ddt(Sxz_)+fvc::div(this->phi(),Sxz_)))/(2.0*magSqr(S));  
    // SDS.component(tensor::YZ) = (Syx_*(fvc::ddt(Szx_)+fvc::div(this->phi(),Szx_))
    //                           + Syy_*(fvc::ddt(Szy_)+fvc::div(this->phi(),Szy_))
    //                           + Syz_*(fvc::ddt(Szz_)+fvc::div(this->phi(),Szz_)))/(2.0*magSqr(S));  
    // SDS.component(tensor::ZX) = (Szx_*(fvc::ddt(Sxx_)+fvc::div(this->phi(),Sxx_))
    //                           + Szy_*(fvc::ddt(Sxy_)+fvc::div(this->phi(),Sxy_))
    //                           + Szz_*(fvc::ddt(Sxz_)+fvc::div(this->phi(),Sxz_)))/(2.0*magSqr(S));  
    // SDS.component(tensor::ZY) = (Szx_*(fvc::ddt(Syx_)+fvc::div(this->phi(),Syx_))
    //                           + Szy_*(fvc::ddt(Syy_)+fvc::div(this->phi(),Syy_))
    //                           + Szz_*(fvc::ddt(Syz_)+fvc::div(this->phi(),Syz_)))/(2.0*magSqr(S));
    
    
    // Spalart-Shur curvature correction for vorticity tensor (TLLP:Eq.20)
    // const volTensorField WTilde
    // (
    //     W - 2.0*skew(SDS)
    // );

    // da testare
    SDS.component(tensor::XY) = (
                                Sxx_*(fvc::ddt(Syx_) + (U & fvc::grad(Syx_)))
                              + Sxy_*(fvc::ddt(Syy_) + (U & fvc::grad(Syy_)))
                              + Sxz_*(fvc::ddt(Syz_) + (U & fvc::grad(Syz_)))
                                )/(2.0*magSqr(S));
    SDS.component(tensor::XZ) = (
                                Sxx_*(fvc::ddt(Szx_) + (U & fvc::grad(Szx_)))
                              + Sxy_*(fvc::ddt(Szy_) + (U & fvc::grad(Szy_)))
                              + Sxz_*(fvc::ddt(Szz_) + (U & fvc::grad(Szz_)))
                                )/(2.0*magSqr(S));                          
    SDS.component(tensor::YX) = (
                                Syx_*(fvc::ddt(Sxx_) + (U & fvc::grad(Sxx_)))
                              + Syy_*(fvc::ddt(Sxy_) + (U & fvc::grad(Sxy_)))
                              + Syz_*(fvc::ddt(Sxz_) + (U & fvc::grad(Sxz_)))
                                )/(2.0*magSqr(S));  
    SDS.component(tensor::YZ) = (
                                Syx_*(fvc::ddt(Szx_) + (U & fvc::grad(Szx_)))
                              + Syy_*(fvc::ddt(Szy_) + (U & fvc::grad(Szy_)))
                              + Syz_*(fvc::ddt(Szz_) + (U & fvc::grad(Szz_)))
                                )/(2.0*magSqr(S));  
    SDS.component(tensor::ZX) = (
                                Szx_*(fvc::ddt(Sxx_) + (U & fvc::grad(Sxx_)))
                              + Szy_*(fvc::ddt(Sxy_) + (U & fvc::grad(Sxy_)))
                              + Szz_*(fvc::ddt(Sxz_) + (U & fvc::grad(Sxz_)))
                                )/(2.0*magSqr(S));  
    SDS.component(tensor::ZY) = (
                                Szx_*(fvc::ddt(Syx_) + (U & fvc::grad(Syx_)))
                              + Szy_*(fvc::ddt(Syy_) + (U & fvc::grad(Syy_)))
                              + Szz_*(fvc::ddt(Syz_) + (U & fvc::grad(Syz_)))
                                )/(2.0*magSqr(S));
    
    // Spalart-Shur curvature correction for vorticity tensor (TLLP:Eq.20)
    const volTensorField WTilde
    (
        W - 2.0*skew(SDS)
    );

    // Anisotropy tensor (TLLP:Eq.18)
    const dimensionedScalar beta2_ = 2.0*(1.0 - C5_)/(C1_ + C1s_ + 1.0);
    volTensorField A
    (
        -2*nut/k_*(S + 2.0*beta2_*((S & WTilde) - (WTilde & S))/
        (mag(S + WTilde)))
    );

    // Time scale
    const volScalarField tau
    (
        k_/epsilon_
    );

    // Function for damping phit in region of low strain (TLLP:Eq.17)
    const volScalarField fmu
    (
        (sqrt(2.0)*mag(S)*tau + pow3(ebf_)) /
        max(
            sqrt(2.0)*mag(S)*tau,1.87
           )
    );

    // Coefficient using fmu (TLLP:Eq.15)
    const volScalarField Cp3
    (
         fmu/Cmu_*(2.0/3.0 - C3_/2.0)
    );

    // Wall-normal vectors defined through the elliptic blending factor
    const volVectorField n 
    (
        fvc::grad(ebf_)/max(
            mag(fvc::grad(ebf_)), dimensionedScalar(dimless/dimLength, SMALL)
            )
    );

    // Additional production term in epsilon eq. (TLLP:Eq.7)
    const volScalarField E
    (
        CK_*pow3(1.0 - ebf_)*this->nu()*nut*sqr(fvc::div(mag(2.0*S & n)*n))
    );
    
    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Turbulent kinetic energy dissipation rate equation  (TLLP:Eq.6)
    // k/T ~ epsilon
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
      ==
      - fvm::SuSp
        (   
            alpha()*rho()*Ceps1_*
            (
                - G()/k_()
                + 2.0/3.0*divU 
            ),
            epsilon_
        )
      - fvm::Sp(alpha()*rho()*Ceps2_/tau(), epsilon_)
      + alpha()*rho()*E()
      + fvOptions(alpha, rho, epsilon_)
    );
    
    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);


    // Turbulent kinetic energy equation (TLLP:Eq.5)
    // epsilon/k ~ 1/Ts
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
      ==
        alpha()*rho()*G()
      - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()/tau(), k_)
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);


// Elliptic blending function equation (TLLP:Eq.11)
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

    // Coefficients to be used in the phitEquation
    const dimensionedScalar Cws = Ceps2_ - 1.0 + 5.0 - 1.0/Cmu_;
    const dimensionedScalar C1Tilde = C1_ + Ceps2_ - 2.0;
    const dimensionedScalar Cp1 = 2.0 - Ceps1_;
    const dimensionedScalar Cp2 = C3s_/sqrt(2.0);
    const dimensionedScalar C4s = 2.0/Cmu_*(1.0 - C4_);
    const dimensionedScalar C5s = 2.0/Cmu_*(1.0 - C5_);

// Stress-strain lag equation (TLLP:Eq.9)
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
                (1.0 - pow3(ebf_()))*Cws/tau()
                + pow3(ebf_())*(C1Tilde + C1s_*(G() - (2.0/3.0)*k_()*divU)
                /epsilon_())/tau()
                + Cp1*(G()/k_() - (2.0/3.0)*divU)
                - pow3(ebf_())*Cp2*sqrt(2.0)*mag(S())
            )
          , phit_
        )
        + alpha()*rho()*
        (
            pow3(ebf_())/tau()/(2.0*magSqr(S()))*((C4s*(A() & S()) 
                - C5s*(A() & WTilde())) && S())
            + pow3(ebf_())*Cp3()/tau()
        )
      + fvOptions(alpha, rho, phit_)
    );

    phitEqn.ref().relax();
    fvOptions.constrain(phitEqn.ref());
    solve(phitEqn);
    fvOptions.correct(phit_);
    bound(phit_, phitMin_);
    
    // Bounding phit
    
    forAll(phit_, celli)
    {
        if(phit_[celli] > 2.0)
        {
            phit_[celli] = 2.0;
        }
    }

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
