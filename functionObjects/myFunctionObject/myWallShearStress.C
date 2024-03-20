/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "myWallShearStress.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(myWallShearStress, 0);
    addToRunTimeSelectionTable(functionObject, myWallShearStress, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::myWallShearStress::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Wall shear stress");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os  << endl;
}


void Foam::functionObjects::myWallShearStress::calcShearStress
(
    const volSymmTensorField& Reff,
    volVectorField& shearStress
)
{
    shearStress.dimensions().reset(Reff.dimensions());

    for (const label patchi : patchSet_)
    {
        vectorField& ssp = shearStress.boundaryFieldRef()[patchi];
        const vectorField& Sfp = mesh_.Sf().boundaryField()[patchi];
        const scalarField& magSfp = mesh_.magSf().boundaryField()[patchi];
        const symmTensorField& Reffp = Reff.boundaryField()[patchi];

        ssp = (-Sfp/magSfp) & Reffp;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::myWallShearStress::myWallShearStress
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name, typeName, dict),
    patchSet_()
{
    read(dict);

    writeFileHeader(file());

    volVectorField* wallShearStressPtr
    (
        new volVectorField
        (
            IOobject
            (
                scopedName(typeName),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(sqr(dimLength)/sqr(dimTime), Zero)
        )
    );

    mesh_.objectRegistry::store(wallShearStressPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::myWallShearStress::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.getOrDefault<wordRes>("patches", wordRes())
        );

    Info<< type() << " " << name() << ":" << nl;

    if (patchSet_.empty())
    {
        forAll(pbm, patchi)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                patchSet_.insert(patchi);
            }
        }

        Info<< "    processing all wall patches" << nl << endl;
    }
    else
    {
        Info<< "    processing wall patches: " << nl;
        labelHashSet filteredPatchSet;
        for (const label patchi : patchSet_)
        {
            if (isA<wallPolyPatch>(pbm[patchi]))
            {
                filteredPatchSet.insert(patchi);
                Info<< "        " << pbm[patchi].name() << endl;
            }
            else
            {
                WarningInFunction
                    << "Requested wall shear stress on non-wall boundary "
                    << "type patch: " << pbm[patchi].name() << endl;
            }
        }

        Info<< endl;

        patchSet_ = filteredPatchSet;
    }

    return true;
}


bool Foam::functionObjects::myWallShearStress::execute()
{
    auto& myWallShearStress =
        mesh_.lookupObjectRef<volVectorField>(scopedName(typeName));

    // Compressible
    {
        typedef compressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcShearStress(modelPtr->devRhoReff(), myWallShearStress);
            return true;
        }
    }

    // Incompressible
    {
        typedef incompressible::turbulenceModel turbType;

        const turbType* modelPtr =
            findObject<turbType>(turbulenceModel::propertiesName);

        if (modelPtr)
        {
            calcShearStress(modelPtr->devReff(), myWallShearStress);
            return true;
        }
    }

    FatalErrorInFunction
        << "Unable to find turbulence model in the "
        << "database" << exit(FatalError);

    return false;
}


bool Foam::functionObjects::myWallShearStress::write()
{
    const auto& myWallShearStress =
        obr_.lookupObject<volVectorField>(scopedName(typeName));

    Log << type() << " " << name() << " write:" << nl
        << "    writing field " << myWallShearStress.name() << endl;

    myWallShearStress.write();

    const fvPatchList& patches = mesh_.boundary();

    for (const label patchi : patchSet_)
    {
        const fvPatch& pp = patches[patchi];

        const vectorField& ssp = myWallShearStress.boundaryField()[patchi];

        vector minSsp = gMin(ssp);
        vector maxSsp = gMax(ssp);
        vector avgSsp = gAverage(ssp);

        if (Pstream::master())
        {
            writeCurrentTime(file());

            file()
                << token::TAB << pp.name()
                << token::TAB << minSsp
                << token::TAB << maxSsp
                << token::TAB << avgSsp
                << endl;
        }

        Log << "    min/max/average(" << pp.name() << ") = "
            << minSsp << ", " << maxSsp  << ", " << avgSsp <<  endl;
    }

    return true;
}


// ************************************************************************* //
