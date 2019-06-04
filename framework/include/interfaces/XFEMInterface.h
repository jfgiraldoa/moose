//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ConsoleStreamInterface.h"
#include "MooseTypes.h"
#include "InputParameters.h"
#include "MooseMesh.h"

class MooseApp;
class AuxiliarySystem;
class NonlinearSystemBase;
class MaterialData;
class FEProblemBase;

namespace libMesh
{
class MeshBase;
class QBase;
} // namespace libMesh

/**
 * This is the XFEMInterface class.  This is an abstract base
 * class that defines interfaces with a class that dynamically
 * modifies the mesh in support of a phantom node approach for XFEM
 */

// ------------------------------------------------------------
// XFEMInterface class definition
class XFEMInterface : public ConsoleStreamInterface
{
public:
  /**
   * Constructor
   */
  explicit XFEMInterface(const InputParameters & params)
    : ConsoleStreamInterface(*params.getCheckedPointerParam<MooseApp *>("_moose_app")),
      _fe_problem(params.getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")),
      _material_data(nullptr),
      _bnd_material_data(nullptr),
      _moose_mesh(nullptr),
      _moose_displaced_mesh(nullptr),
      _mesh(nullptr),
      _displaced_mesh(nullptr)
  {
  }

  /**
   * Destructor
   */
  virtual ~XFEMInterface() {}

  /**
   * Set the pointer to the primary mesh that is modified by XFEM
   */
  void setMesh(MooseMesh * mesh)
  {
    _moose_mesh = mesh;
    _mesh = &mesh->getMesh();
  }

  /**
   * Set the pointer to the displaced mesh that is modified by XFEM
   */
  void setDisplacedMesh(MooseMesh * displaced_mesh)
  {
    _moose_displaced_mesh = displaced_mesh;
    _displaced_mesh = &displaced_mesh->getMesh();
  }

  /**
   * Set the pointer to the MaterialData
   */
  void setMaterialData(std::vector<std::shared_ptr<MaterialData>> * material_data)
  {
    _material_data = material_data;
  }

  /**
   * Set the pointer to the Boundary MaterialData
   */
  void setBoundaryMaterialData(std::vector<std::shared_ptr<MaterialData>> * bnd_material_data)
  {
    _bnd_material_data = bnd_material_data;
  }

  /**
   * Method to update the mesh due to modified cut definitions
   */
  virtual bool update(Real time, NonlinearSystemBase & nl, AuxiliarySystem & aux) = 0;

  /**
   * Initialize the solution on newly created nodes
   */
  virtual void initSolution(NonlinearSystemBase & nl, AuxiliarySystem & aux) = 0;

  /**
   * Get the factors for the QP weighs for XFEM partial elements
   * @param weights The new weights at element quadrature points
   * @param elem The element for which the weights are adjusted
   * @param qrule The quadrature rule for the volume integration
   * @param q_points The vector of quadrature points
   */

  virtual bool getXFEMWeights(MooseArray<Real> & weights,
                              const Elem * elem,
                              QBase * qrule,
                              const MooseArray<Point> & q_points) = 0;

  /**
   * Get the factors for the face QP weighs for XFEM partial elements
   * @param weights The new weights at element face quadrature points
   * @param elem The element for which the weights are adjusted
   * @param qrule The quadrature rule for the face integration
   * @param q_points The vector of quadrature points at element face
   * @param side The side of element for which the weights are adjusted
   */
  virtual bool getXFEMFaceWeights(MooseArray<Real> & weights,
                                  const Elem * elem,
                                  QBase * qrule,
                                  const MooseArray<Point> & q_points,
                                  unsigned int side) = 0;

  /**
   * Potentially update the mesh by healing previous XFEM cuts.
   * @return true if the mesh has been updated due to healing
   **/
  virtual bool updateHeal() = 0;

protected:
  FEProblemBase * _fe_problem;
  std::vector<std::shared_ptr<MaterialData>> * _material_data;
  std::vector<std::shared_ptr<MaterialData>> * _bnd_material_data;

  MooseMesh * _moose_mesh;
  MooseMesh * _moose_displaced_mesh;
  MeshBase * _mesh;
  MeshBase * _displaced_mesh;
};

