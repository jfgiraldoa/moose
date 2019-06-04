//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "MultiAppFieldTransfer.h"

// Forward declarations
class MultiAppNearestNodeTransfer;
namespace libMesh
{
class DofObject;
}

template <>
InputParameters validParams<MultiAppNearestNodeTransfer>();

/**
 * Copy the value to the target domain from the nearest node in the source domain.
 */
class MultiAppNearestNodeTransfer : public MultiAppFieldTransfer
{
public:
  MultiAppNearestNodeTransfer(const InputParameters & parameters);

  virtual void execute() override;

protected:
  /**
   * Return the nearest node to the point p.
   * @param p The point you want to find the nearest node to.
   * @param distance This will hold the distance between the returned node and p
   * @param mesh The mesh in which we search for the node
   * @param local true if we look at local nodes, otherwise we look at all nodes
   * @return The Node closest to point p.
   */
  Node * getNearestNode(const Point & p, Real & distance, MooseMesh * mesh, bool local);

  /**
   * Return the distance between the given point and the farthest corner of the
   * given bounding box.
   * @param p The point to evaluate all distances from.
   * @param bbox The bounding box to evaluate the distance to.
   * @return The maximum distance between the point p and the eight corners of
   * the bounding box bbox.
   */
  Real bboxMaxDistance(const Point & p, const BoundingBox & bbox);

  /**
   * Return the distance between the given point and the nearest corner of the
   * given bounding box.
   * @param p The point to evaluate all distances from.
   * @param bbox The bounding box to evaluate the distance to.
   * @return The minimum distance between the point p and the eight corners of
   * the bounding box bbox.
   */
  Real bboxMinDistance(const Point & p, const BoundingBox & bbox);

  void getLocalEntities(MooseMesh * mesh,
                        std::vector<std::pair<Point, DofObject *>> & local_entities,
                        bool nodal);

  /// If true then node connections will be cached
  bool _fixed_meshes;

  /// Used to cache nodes
  std::map<dof_id_type, Node *> & _node_map;

  /// Used to cache distances
  std::map<dof_id_type, Real> & _distance_map;

  // These variables allow us to cache nearest node info
  bool & _neighbors_cached;
  std::vector<std::vector<unsigned int>> & _cached_froms;
  std::vector<std::vector<dof_id_type>> & _cached_dof_ids;
  std::map<dof_id_type, unsigned int> & _cached_from_inds;
  std::map<dof_id_type, unsigned int> & _cached_qp_inds;
};

