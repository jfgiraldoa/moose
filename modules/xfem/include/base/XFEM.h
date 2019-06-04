//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementPairLocator.h"
#include "ElementFragmentAlgorithm.h"
#include "XFEMInterface.h"
#include "XFEMCrackGrowthIncrement2DCut.h"

#include "libmesh/vector_value.h"
#include "libmesh/quadrature.h"

#include "GeometricCutUserObject.h"

// Forward declarations
class SystemBase;

namespace Xfem
{
enum XFEM_CUTPLANE_QUANTITY
{
  ORIGIN_X,
  ORIGIN_Y,
  ORIGIN_Z,
  NORMAL_X,
  NORMAL_Y,
  NORMAL_Z
};

enum XFEM_QRULE
{
  VOLFRAC,
  MOMENT_FITTING,
  DIRECT
};
} // namespace Xfem

class XFEMCutElem;
class XFEMCrackGrowthIncrement2DCut;
class EFANode;
class EFAEdge;
class EFAElement;
class EFAElement2D;
class EFAElement3D;

/**
 * This is the \p XFEM class.  This class implements
 * algorithms for dynamic mesh modification in support of
 * a phantom node approach for XFEM
 */

// ------------------------------------------------------------
// XFEM class definition
class XFEM : public XFEMInterface
{
public:
  explicit XFEM(const InputParameters & params);

  ~XFEM();

  void addGeometricCut(GeometricCutUserObject * geometric_cut);

  void addStateMarkedElem(unsigned int elem_id, RealVectorValue & normal);
  void addStateMarkedElem(unsigned int elem_id, RealVectorValue & normal, unsigned int marked_side);
  void addStateMarkedFrag(unsigned int elem_id, RealVectorValue & normal);

  void clearStateMarkedElems();

  /**
   * Add information about a new cut to be performed on a specific 2d element
   * @param elem_id   The id of the element to be cut
   * @param geom_info The object containing information about the cut to be performed
   * @param interface_id The ID of the interface
   */
  void addGeomMarkedElem2D(const unsigned int elem_id,
                           const Xfem::GeomMarkedElemInfo2D geom_info,
                           const unsigned int interface_id);

  /**
   * Add information about a new cut to be performed on a specific 3d element
   * @param elem_id   The id of the element to be cut
   * @param geom_info The object containing information about the cut to be performed
   * @param interface_id The ID of the interface
   */
  void addGeomMarkedElem3D(const unsigned int elem_id,
                           const Xfem::GeomMarkedElemInfo3D geom_info,
                           const unsigned int interface_id);

  /**
   * Clear out the list of elements to be marked for cutting. Called after cutting is done.
   */
  void clearGeomMarkedElems();

  virtual bool update(Real time, NonlinearSystemBase & nl, AuxiliarySystem & aux) override;

  virtual void initSolution(NonlinearSystemBase & nl, AuxiliarySystem & aux) override;

  void buildEFAMesh();
  bool markCuts(Real time);
  bool markCutEdgesByGeometry();
  bool markCutEdgesByState(Real time);
  bool markCutFacesByGeometry();
  bool markCutFacesByState();
  bool initCutIntersectionEdge(
      Point cut_origin, RealVectorValue cut_normal, Point & edge_p1, Point & edge_p2, Real & dist);
  bool cutMeshWithEFA(NonlinearSystemBase & nl, AuxiliarySystem & aux);

  /**
   * Potentially heal the mesh by merging some of the pairs
   * of partial elements cut by XFEM back into single elements
   * if indicated by the cutting objects.
   * @return true if the mesh has been modified due to healing
   **/
  bool healMesh();

  virtual bool updateHeal() override;
  Point getEFANodeCoords(EFANode * CEMnode,
                         EFAElement * CEMElem,
                         const Elem * elem,
                         MeshBase * displaced_mesh = NULL) const;

  /**
   * Get the volume fraction of an element that is physical
   */
  Real getPhysicalVolumeFraction(const Elem * elem) const;

  /**
   * Return true if the point is inside the element physical domain
   * Note: if this element is not cut, return true too
   */
  bool isPointInsidePhysicalDomain(const Elem * elem, const Point & point) const;

  /**
   * Get specified component of normal or origin for cut plane for a given element
   */
  Real getCutPlane(const Elem * elem,
                   const Xfem::XFEM_CUTPLANE_QUANTITY quantity,
                   unsigned int plane_id) const;

  bool isElemAtCrackTip(const Elem * elem) const;
  bool isElemCut(const Elem * elem, XFEMCutElem *& xfce) const;
  bool isElemCut(const Elem * elem) const;
  void getFragmentFaces(const Elem * elem,
                        std::vector<std::vector<Point>> & frag_faces,
                        bool displaced_mesh = false) const;
  void storeCrackTipOriginAndDirection();

  void correctCrackExtensionDirection(const Elem * elem,
                                      EFAElement2D * CEMElem,
                                      EFAEdge * orig_edge,
                                      Point normal,
                                      Point crack_tip_origin,
                                      Point crack_tip_direction,
                                      Real & distance_keep,
                                      unsigned int & edge_id_keep,
                                      Point & normal_keep);

  void getCrackTipOrigin(std::map<unsigned int, const Elem *> & elem_id_crack_tip,
                         std::vector<Point> & crack_front_points);

  Xfem::XFEM_QRULE & getXFEMQRule();
  void setXFEMQRule(std::string & xfem_qrule);
  void setCrackGrowthMethod(bool use_crack_growth_increment, Real crack_growth_increment);
  virtual bool getXFEMWeights(MooseArray<Real> & weights,
                              const Elem * elem,
                              QBase * qrule,
                              const MooseArray<Point> & q_points) override;
  virtual bool getXFEMFaceWeights(MooseArray<Real> & weights,
                                  const Elem * elem,
                                  QBase * qrule,
                                  const MooseArray<Point> & q_points,
                                  unsigned int side) override;

  /**
   * Get the list of cut element pairs corresponding to a given
   * interface ID.
   * @param interface_id The ID of the interface
   * @return the list of elements cut by that interface
   **/
  virtual const ElementPairLocator::ElementPairList * getXFEMCutElemPairs(unsigned int interface_id)
  {
    return &_sibling_elems[interface_id];
  }

  /**
   * Get the list of cut element pairs on the displaced mesh
   * corresponding to a given interface ID.
   * @param interface_id The ID of the interface
   * @return the list of elements cut by that interface
   **/
  virtual const ElementPairLocator::ElementPairList *
  getXFEMDisplacedCutElemPairs(unsigned int interface_id)
  {
    return &_sibling_displaced_elems[interface_id];
  }

  /**
   * Get the interface ID corresponding to a given GeometricCutUserObject.
   * @param gcu pointer to the GeometricCutUserObject
   * @return the interface ID
   **/
  virtual unsigned int getGeometricCutID(const GeometricCutUserObject * gcu)
  {
    return _geom_marker_id_map[gcu];
  };

  virtual void getXFEMIntersectionInfo(const Elem * elem,
                                       unsigned int plane_id,
                                       Point & normal,
                                       std::vector<Point> & intersectionPoints,
                                       bool displaced_mesh = false) const;
  virtual void getXFEMqRuleOnLine(std::vector<Point> & intersection_points,
                                  std::vector<Point> & quad_pts,
                                  std::vector<Real> & quad_wts) const;
  virtual void getXFEMqRuleOnSurface(std::vector<Point> & intersection_points,
                                     std::vector<Point> & quad_pts,
                                     std::vector<Real> & quad_wts) const;
  bool has_secondary_cut() { return _has_secondary_cut; }

  /**
   * Get the EFAElement2D object for a specified libMesh element
   * @param elem Pointer to the libMesh element for which the object is requested
   */
  EFAElement2D * getEFAElem2D(const Elem * elem);

  /**
   * Get the EFAElement3D object for a specified libMesh element
   * @param elem Pointer to the libMesh element for which the object is requested
   */
  EFAElement3D * getEFAElem3D(const Elem * elem);

  void getFragmentEdges(const Elem * elem,
                        EFAElement2D * CEMElem,
                        std::vector<std::vector<Point>> & frag_edges) const;
  void getFragmentFaces(const Elem * elem,
                        EFAElement3D * CEMElem,
                        std::vector<std::vector<Point>> & frag_faces) const;

  const std::map<const Elem *, std::vector<Point>> & getCrackTipOriginMap() const
  {
    return _elem_crack_origin_direction_map;
  }

private:
  bool _has_secondary_cut;

  Xfem::XFEM_QRULE _XFEM_qrule;

  bool _use_crack_growth_increment;
  Real _crack_growth_increment;

  std::vector<const GeometricCutUserObject *> _geometric_cuts;

  std::map<unique_id_type, XFEMCutElem *> _cut_elem_map;
  std::set<const Elem *> _crack_tip_elems;
  std::set<const Elem *> _crack_tip_elems_to_be_healed;
  std::map<unsigned int, ElementPairLocator::ElementPairList> _sibling_elems;
  std::map<unsigned int, ElementPairLocator::ElementPairList> _sibling_displaced_elems;

  std::map<const Elem *, std::vector<Point>> _elem_crack_origin_direction_map;

  // std::map<const Elem*, Point> _crack_propagation_direction_map;

  std::map<const Elem *, RealVectorValue> _state_marked_elems;
  std::set<const Elem *> _state_marked_frags;
  std::map<const Elem *, unsigned int> _state_marked_elem_sides;

  /// Data structure for storing information about all 2D elements to be cut by geometry
  std::map<const Elem *, std::vector<Xfem::GeomMarkedElemInfo2D>> _geom_marked_elems_2d;

  /// Data structure for storing information about all 3D elements to be cut by geometry
  std::map<const Elem *, std::vector<Xfem::GeomMarkedElemInfo3D>> _geom_marked_elems_3d;

  /// Data structure for storing the elements cut by specific geometric cutters
  std::map<unsigned int, std::set<unsigned int>> _geom_marker_id_elems;

  /// Data structure for storing the GeommetricCutUserObjects and their corresponding id
  std::map<const GeometricCutUserObject *, unsigned int> _geom_marker_id_map;

  ElementFragmentAlgorithm _efa_mesh;

  /**
   * Data structure to store the nonlinear solution for nodes/elements affected by XFEM
   * For each node/element, this is stored as a vector that contains all components
   * of all applicable variables in an order defined by getElementSolutionDofs() or
   * getNodeSolutionDofs(). This vector first contains the current solution in that
   * order, followed by the old and older solutions, also in that same order.
   */
  std::map<unique_id_type, std::vector<Real>> _cached_solution;

  /**
   * Data structure to store the auxiliary solution for nodes/elements affected by XFEM
   * For each node/element, this is stored as a vector that contains all components
   * of all applicable variables in an order defined by getElementSolutionDofs() or
   * getNodeSolutionDofs(). This vector first contains the current solution in that
   * order, followed by the old and older solutions, also in that same order.
   */
  std::map<unique_id_type, std::vector<Real>> _cached_aux_solution;

  /**
   * Store the solution in stored_solution for a given node
   * @param node_to_store_to   Node for which the solution will be stored
   * @param node_to_store_from Node from which the solution to be stored is obtained
   * @param sys                System from which the solution is stored
   * @param stored_solution    Data structure that the stored solution is saved to
   * @param current_solution   Current solution vector that the solution is obtained from
   * @param old_solution       Old solution vector that the solution is obtained from
   * @param older_solution     Older solution vector that the solution is obtained from
   */
  void storeSolutionForNode(const Node * node_to_store_to,
                            const Node * node_to_store_from,
                            SystemBase & sys,
                            std::map<unique_id_type, std::vector<Real>> & stored_solution,
                            const NumericVector<Number> & current_solution,
                            const NumericVector<Number> & old_solution,
                            const NumericVector<Number> & older_solution);

  /**
   * Store the solution in stored_solution for a given element
   * @param elem_to_store_to   Element for which the solution will be stored
   * @param elem_to_store_from Element from which the solution to be stored is obtained
   * @param sys                System from which the solution is stored
   * @param stored_solution    Data structure that the stored solution is saved to
   * @param current_solution   Current solution vector that the solution is obtained from
   * @param old_solution       Old solution vector that the solution is obtained from
   * @param older_solution     Older solution vector that the solution is obtained from
   */
  void storeSolutionForElement(const Elem * elem_to_store_to,
                               const Elem * elem_to_store_from,
                               SystemBase & sys,
                               std::map<unique_id_type, std::vector<Real>> & stored_solution,
                               const NumericVector<Number> & current_solution,
                               const NumericVector<Number> & old_solution,
                               const NumericVector<Number> & older_solution);

  /**
   * Set the solution for all locally-owned nodes/elements that have stored values
   * @param sys              System for which the solution is set
   * @param stored_solution  Data structure that the stored solution is obtained from
   * @param current_solution Current solution vector that will be set
   * @param old_solution     Old solution vector that will be set
   * @param older_solution   Older solution vector that will be set
   */
  void setSolution(SystemBase & sys,
                   const std::map<unique_id_type, std::vector<Real>> & stored_solution,
                   NumericVector<Number> & current_solution,
                   NumericVector<Number> & old_solution,
                   NumericVector<Number> & older_solution);

  /**
   * Set the solution for a set of DOFs
   * @param stored_solution      Stored solution values to set the solution to
   * @param stored_solution_dofs Dof indices for the entries in stored_solution
   * @param current_solution     Current solution vector that will be set
   * @param old_solution         Old solution vector that will be set
   * @param older_solution       Older solution vector that will be set
   */
  void setSolutionForDOFs(const std::vector<Real> & stored_solution,
                          const std::vector<dof_id_type> & stored_solution_dofs,
                          NumericVector<Number> & current_solution,
                          NumericVector<Number> & old_solution,
                          NumericVector<Number> & older_solution);

  /**
   * Get a vector of the dof indices for all components of all variables
   * associated with an element
   * @param elem Element for which dof indices are found
   * @param sys  System for which the dof indices are found
   */
  std::vector<dof_id_type> getElementSolutionDofs(const Elem * elem, SystemBase & sys) const;

  /**
   * Get a vector of the dof indices for all components of all variables
   * associated with a node
   * @param node Node for which dof indices are found
   * @param sys  System for which the dof indices are found
   */
  std::vector<dof_id_type> getNodeSolutionDofs(const Node * node, SystemBase & sys) const;
};

