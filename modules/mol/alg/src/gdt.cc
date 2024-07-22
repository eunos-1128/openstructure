//------------------------------------------------------------------------------
// This file is part of the OpenStructure project <www.openstructure.org>
//
// Copyright (C) 2008-2024 by the OpenStructure authors
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 3.0 of the License, or (at your option)
// any later version.
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
//------------------------------------------------------------------------------

#include <ost/mol/alg/gdt.hh>
#include <ost/message.hh>
#include <Eigen/Dense>

namespace {

// get RMSD and rotation matrix using the Theobald method
// Both position matrices are expected to have the same size
// and to have their average position at the origin 
void TheobaldRMSD(const Eigen::Matrix<double,Eigen::Dynamic,3>& pos_one,
                  const Eigen::Matrix<double,Eigen::Dynamic,3>& pos_two,
                  Real& rmsd, Eigen::Matrix<double,3,3>& rot){

  if(pos_one.rows() < 3){
    throw ost::Error("Observed superposition with < 3 positions to "
                         "superpose!");
  }

  Eigen::Matrix<double,3,3> M = pos_one.transpose() * pos_two;

  // using floats for the squared norm is fine
  double GA = pos_one.squaredNorm();
  double GB = pos_two.squaredNorm();

  Eigen::Matrix<double,4,4> K;
  K(0,0)  =  M(0,0) + M(1,1) + M(2,2);
  K(0,1)  =  M(1,2) - M(2,1);
  K(0,2)  =  M(2,0) - M(0,2);
  K(0,3)  =  M(0,1) - M(1,0);

  K(1,0)  =  K(0,1);
  K(1,1)  =  M(0,0) - M(1,1) - M(2,2);
  K(1,2)  =  M(0,1) + M(1,0);
  K(1,3)  =  M(0,2) + M(2,0);

  K(2,0)  =  K(0,2);
  K(2,1)  =  K(1,2);
  K(2,2)  = -M(0,0) + M(1,1) - M(2,2);
  K(2,3)  =  M(1,2) + M(2,1);

  K(3,0)  =  K(0,3);
  K(3,1)  =  K(1,3);
  K(3,2)  =  K(2,3);
  K(3,3)  = -M(0,0) - M(1,1) + M(2,2);

  double C0 = K.determinant();
  double C1 = -8.0*M.determinant();
  double C2 = -2.0*M.squaredNorm();
  double lambda = 0.5 * (GA + GB);
  double a, b, d, lambda_2;
  for(int i = 0; i < 50; ++i){
    lambda_2 = lambda * lambda;
    b = (lambda_2 + C2) * lambda;
    a = b + C1;
    d = (a*lambda + C0) / (2.0*lambda_2*lambda + b + a);
    lambda -= d;
    if(std::abs(d) < 1e-6){
      break;
    }
  }

  double msd = (GA + GB - 2.0 * lambda) / pos_one.rows();
  if(msd < 1e-4){
    // The algorithm never really goes to zero... if msd is super small we just
    // assign zero. 1e-4 corresponds to an rmsd of 0.01
    rmsd = 0.0;
  }
  else{
    rmsd = std::sqrt(msd);
  }

  K -= lambda*Eigen::Matrix<double,4,4>::Identity();

  double helper[6];
  helper[0] = K(2,2)*K(3,3) - K(3,2)*K(2,3);
  helper[1] = K(2,1)*K(3,3) - K(3,1)*K(2,3);
  helper[2] = K(2,1)*K(3,2) - K(3,1)*K(2,2);
  helper[3] = K(2,0)*K(3,3) - K(3,0)*K(2,3);
  helper[4] = K(2,0)*K(3,2) - K(3,0)*K(2,2);
  helper[5] = K(2,0)*K(3,1) - K(3,0)*K(2,1);

  double q1 =  K(1,1)*helper[0] - K(1,2)*helper[1] + K(1,3)*helper[2];
  double q2 = -K(1,0)*helper[0] + K(1,2)*helper[3] - K(1,3)*helper[4];
  double q3 =  K(1,0)*helper[1] - K(1,1)*helper[3] + K(1,3)*helper[5];
  double q4 = -K(1,0)*helper[2] + K(1,1)*helper[4] - K(1,2)*helper[5]; 
  double norm = q1*q1 + q2*q2 + q3*q3 + q4*q4;

  if(norm < 1e-6){
    q1 =  K(0,1)*helper[0] - K(0,2)*helper[1] + K(0,3)*helper[2];
    q2 = -K(0,0)*helper[0] + K(0,2)*helper[3] - K(0,3)*helper[4];
    q3 =  K(0,0)*helper[1] - K(0,1)*helper[3] + K(0,3)*helper[5];
    q4 = -K(0,0)*helper[2] + K(0,1)*helper[4] - K(0,2)*helper[5];
    norm = q1*q1 + q2*q2 +q3*q3 + q4*q4;

    if (norm < 1e-6){
      helper[0] = K(0,2)*K(1,3) - K(0,3)*K(1,2);
      helper[1] = K(0,1)*K(1,3) - K(0,3)*K(1,1);
      helper[2] = K(0,1)*K(1,2) - K(0,2)*K(1,1);
      helper[3] = K(0,0)*K(1,3) - K(0,3)*K(1,0);
      helper[4] = K(0,0)*K(1,2) - K(0,2)*K(1,0);
      helper[5] = K(0,0)*K(1,1) - K(0,1)*K(1,0);

      q1 =  K(3,1)*helper[0] - K(3,2)*helper[1] + K(3,3)*helper[2];
      q2 = -K(3,0)*helper[0] + K(3,2)*helper[3] - K(3,3)*helper[4];
      q3 =  K(3,0)*helper[1] - K(3,1)*helper[3] + K(3,3)*helper[5];
      q4 = -K(3,0)*helper[2] + K(3,1)*helper[4] - K(3,2)*helper[5];
      norm = q1*q1 + q2*q2 + q3*q3 + q4*q4;

      if (norm < 1e-6){
        q1 =  K(2,1)*helper[0] - K(2,2)*helper[1] + K(2,3)*helper[2];
        q2 = -K(2,0)*helper[0] + K(2,2)*helper[3] - K(2,3)*helper[4];
        q3 =  K(2,0)*helper[1] - K(2,1)*helper[3] + K(2,3)*helper[5];
        q4 = -K(2,0)*helper[2] + K(2,1)*helper[4] - K(2,2)*helper[5];
        norm = q1*q1 + q2*q2 + q3*q3 + q4*q4;
        if (norm < 1e-6){
          // this should not happen
          rot = Eigen::Matrix<double,3,3>::Identity();
          return;
        }
      }
    }
  }

  norm = 1.0 / std::sqrt(norm);
  q1 *= norm; q2 *= norm; q3 *= norm; q4 *= norm;
  rot(0,0) = 1.0 - 2.0*(q3*q3 + q4*q4);
  rot(0,1) =       2.0*(q2*q3 - q1*q4);
  rot(0,2) =       2.0*(q2*q4 + q3*q1);
  rot(1,0) =       2.0*(q2*q3 + q1*q4);
  rot(1,1) = 1.0 - 2.0*(q2*q2 + q4*q4);
  rot(1,2) =       2.0*(q3*q4 - q2*q1);
  rot(2,0) =       2.0*(q2*q4 - q3*q1);
  rot(2,1) =       2.0*(q3*q4 + q2*q1);
  rot(2,2) = 1.0 - 2.0*(q2*q2 + q3*q3);
}

void Superpose(Eigen::Matrix<double,Eigen::Dynamic,3>& pos_one,
               Eigen::Matrix<double,Eigen::Dynamic,3>& pos_two,
               Eigen::Matrix<double,1,3>& avg_one,
               Eigen::Matrix<double,1,3>& avg_two,
               Real& rmsd,
               Eigen::Matrix<double,3,3>& rotation){

  if(pos_one.rows() != pos_two.rows()){
    throw ost::Error("Cannot superpose positions of different size!");
  }


  avg_one = Eigen::Matrix<double,1,3>::Zero();
  for (uint i = 0; i < pos_one.rows(); ++i) {
    avg_one += pos_one.row(i);
  }
  avg_one = avg_one / pos_one.rows();

  avg_two = Eigen::Matrix<double,1,3>::Zero();
  for (uint i = 0; i < pos_two.rows(); ++i) {
    avg_two += pos_two.row(i);
  }
  avg_two = avg_two / pos_two.rows();

  // TheobaldRMSD only determines the rotational component of the superposition
  // we need to shift the centers of the two point sets onto origin
  for (uint i = 0; i < pos_one.rows(); ++i){
    pos_one.row(i) -= avg_one;
    pos_two.row(i) -= avg_two;
  }

  TheobaldRMSD(pos_one, pos_two, rmsd, rotation);
}

void SuperposeIterative(const Eigen::Matrix<double,Eigen::Dynamic,3>& pos_one,
                        const Eigen::Matrix<double,Eigen::Dynamic,3>& pos_two,
                        int max_iterations, 
                        Real distance_thresh,
                        std::vector<int>& indices,
                        Eigen::Matrix<double,1,3>& avg_one,
                        Eigen::Matrix<double,1,3>& avg_two,
                        Eigen::Matrix<double,3,3>& rotation){

  if(pos_one.rows() != pos_two.rows()){
    throw ost::Error("Position data must be of consistent size!");
  }

  if(max_iterations <= 0){
    throw ost::Error("max_iterations must be at least 1!");
  }

  int num_rows = pos_one.rows();

  if(indices.empty()){
    //there are no idx, so we use all positions for the initial superposition
    indices.resize(num_rows);
    for(int i = 0; i < num_rows; ++i){
      indices[i] = i;
    }
  }
  else{
    //its not empty! let's quickly check whether there are at least 3 positions
    //and whether the indices are all valid
    if(indices.size() < 3){
      throw ost::Error("Must have at least 3 start indices for iterative " 
                       "Superposition!");
    }
    for(auto i: indices) {
      if(i >= num_rows){
        throw ost::Error("Invalid index in iterative Superposition!");
      }
    }
  }

  Real squared_dist_thresh = distance_thresh * distance_thresh;
  std::vector<int> new_indices;
  new_indices.reserve(num_rows);
  std::vector<Real> squared_distances(num_rows);
  Eigen::Matrix<double,1,3> temp_vec;
  Real rmsd;

  // keep track of the indices which give the superposition with maximum
  // number of superposed positions. Thats not necessarily the last superposition
  // when the algorithm converges
  int max_n = -1;
  std::vector<int> temp_indices = indices;
  Eigen::Matrix<double,1,3> temp_avg_one = Eigen::Matrix<double,1,3>::Zero();
  Eigen::Matrix<double,1,3> temp_avg_two = Eigen::Matrix<double,1,3>::Zero();
  Eigen::Matrix<double,3,3> temp_rotation = Eigen::Matrix<double,3,3>::Identity();

  for(int iteration = 0; iteration < max_iterations; ++iteration){

    if(temp_indices.size() < 3) break; //the thing is not really superposable...

    Eigen::Matrix<double,Eigen::Dynamic,3> temp_pos_one = 
    Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(temp_indices.size(), 3);
    Eigen::Matrix<double,Eigen::Dynamic,3> temp_pos_two = 
    Eigen::Matrix<double,Eigen::Dynamic,3>::Zero(temp_indices.size(), 3);

    for(uint i = 0; i < temp_indices.size(); ++i){
      temp_pos_one.row(i) = pos_one.row(temp_indices[i]);
      temp_pos_two.row(i) = pos_two.row(temp_indices[i]);
    }

    Superpose(temp_pos_one, temp_pos_two, temp_avg_one, temp_avg_two, rmsd,
              temp_rotation);

    for(int i = 0; i < num_rows; ++i){
      temp_vec = pos_one.row(i) - temp_avg_one;
      temp_vec = (temp_rotation * temp_vec.transpose()).transpose() + temp_avg_two;
      squared_distances[i] = (temp_vec - pos_two.row(i)).squaredNorm();
    }

    new_indices.clear();
    for(int i = 0; i < num_rows; ++i){
      if(squared_distances[i] < squared_dist_thresh){
        new_indices.push_back(i);
      }
    }

    if(static_cast<int>(new_indices.size()) > max_n) {
      max_n = new_indices.size();
      indices = new_indices;
      avg_one = temp_avg_one;
      avg_two = temp_avg_two;
      rotation = temp_rotation;
    }

    if(new_indices == temp_indices) break; //nothing changes anymore

    temp_indices = new_indices;
  }
}

}

namespace ost{ namespace mol{ namespace alg {


void GDT(const geom::Vec3List& mdl_pos, const geom::Vec3List& ref_pos,
         int window_size, int max_windows, Real distance_thresh,
         int& n_superposed, geom::Mat4& transform) {

  if(mdl_pos.size() != ref_pos.size()){
    throw ost::Error("Position data must be of consistent size!");
  }

  int n_pos = mdl_pos.size();

  // deal with special cases that don't produce valid transforms
  if(n_pos == 1) {
      transform = geom::Mat4::Identity();
      transform.PasteTranslation(ref_pos[0] - mdl_pos[0]);
      n_superposed = 1;
    return;
  }

  if(n_pos == 2) {
    Real mdl_d = geom::Distance(mdl_pos[0], mdl_pos[1]);
    Real ref_d = geom::Distance(ref_pos[0], ref_pos[1]);
    Real dd = std::abs(mdl_d - ref_d);
    if(dd/2 <= distance_thresh) {
      // the two can be superposed within specified distance threshold
      // BUT: cannot construct valid transformation from two positions
      // => Construct matrix with four positions 
      // Two are constructed starting from the center point +- some direction
      // vector that is orthogonal to the vector connecting the original two
      // points. 
      geom::Vec3 mdl_center = (mdl_pos[0] + mdl_pos[1])*0.5;
      geom::Vec3 ref_center = (ref_pos[0] + ref_pos[1])*0.5;
      Eigen::Matrix<double, 4, 3> eigen_mdl_pos = \
      Eigen::Matrix<double, 4, 3>::Zero(4, 3);
      Eigen::Matrix<double, 4, 3> eigen_ref_pos = \
      Eigen::Matrix<double, 4, 3>::Zero(4, 3);
      Eigen::Matrix<double,3,3> eigen_rot = \
      Eigen::Matrix<double,3,3>::Identity();

      geom::Vec3 mdl_dir = geom::Normalize(mdl_pos[1] - mdl_pos[0]);
      geom::Vec3 ref_dir = geom::Normalize(ref_pos[1] - ref_pos[0]);
      geom::Vec3 mdl_normal;
      geom::Vec3 ref_normal;

      // Use cross product to get some normal on mdl_dir
      // The direction of the second vector doesn't really matter, but shouldnt
      // be collinear with mdl_dir
      if(mdl_dir[0] < 0.999) {
        mdl_normal = geom::Cross(geom::Vec3(1,0,0), mdl_dir);
      } else {
        mdl_normal = geom::Cross(geom::Vec3(0,1,0), mdl_dir);
      }

      // same for ref_dir
      if(ref_dir[0] < 0.999) {
        ref_normal = geom::Cross(geom::Vec3(1,0,0), ref_dir);
      } else {
        ref_normal = geom::Cross(geom::Vec3(0,1,0), ref_dir);
      }

      eigen_mdl_pos(0, 0) = mdl_pos[0][0] - mdl_center[0];
      eigen_mdl_pos(0, 1) = mdl_pos[0][1] - mdl_center[1];
      eigen_mdl_pos(0, 2) = mdl_pos[0][2] - mdl_center[2];
      eigen_mdl_pos(1, 0) = mdl_pos[1][0] - mdl_center[0];
      eigen_mdl_pos(1, 1) = mdl_pos[1][1] - mdl_center[1];
      eigen_mdl_pos(1, 2) = mdl_pos[1][2] - mdl_center[2];
      eigen_mdl_pos(2, 0) = mdl_normal[0];
      eigen_mdl_pos(2, 1) = mdl_normal[1];
      eigen_mdl_pos(2, 2) = mdl_normal[2];
      eigen_mdl_pos(3, 0) = -mdl_normal[0];
      eigen_mdl_pos(3, 1) = -mdl_normal[1];
      eigen_mdl_pos(3, 2) = -mdl_normal[2];
      eigen_ref_pos(0, 0) = ref_pos[0][0] - ref_center[0];
      eigen_ref_pos(0, 1) = ref_pos[0][1] - ref_center[1];
      eigen_ref_pos(0, 2) = ref_pos[0][2] - ref_center[2];
      eigen_ref_pos(1, 0) = ref_pos[1][0] - ref_center[0];
      eigen_ref_pos(1, 1) = ref_pos[1][1] - ref_center[1];
      eigen_ref_pos(1, 2) = ref_pos[1][2] - ref_center[2];
      eigen_ref_pos(2, 0) = ref_normal[0];
      eigen_ref_pos(2, 1) = ref_normal[1];
      eigen_ref_pos(2, 2) = ref_normal[2];
      eigen_ref_pos(3, 0) = -ref_normal[0];
      eigen_ref_pos(3, 1) = -ref_normal[1];
      eigen_ref_pos(3, 2) = -ref_normal[2];

      Real tmp; // no need to store RMSD
      TheobaldRMSD(eigen_mdl_pos, eigen_ref_pos, tmp, eigen_rot);
      transform = geom::Mat4();
      transform(0, 0) = eigen_rot(0, 0);
      transform(0, 1) = eigen_rot(0, 1);
      transform(0, 2) = eigen_rot(0, 2);
      transform(1, 0) = eigen_rot(1, 0);
      transform(1, 1) = eigen_rot(1, 1);
      transform(1, 2) = eigen_rot(1, 2);
      transform(2, 0) = eigen_rot(2, 0);
      transform(2, 1) = eigen_rot(2, 1);
      transform(2, 2) = eigen_rot(2, 2);

      // there are three transformation to be applied to reach ref_pos from
      // mdl_pos:
      // 1: shift mdl_pos to center
      // 2: apply estimated rotation
      // 3: shift onto average of ref_pos
      Eigen::Matrix<double,1,3> eigen_avg_mdl = Eigen::Matrix<double,1,3>::Zero();
      Eigen::Matrix<double,1,3> eigen_avg_ref = Eigen::Matrix<double,1,3>::Zero();
      eigen_avg_mdl(0,0) = mdl_center[0]; 
      eigen_avg_mdl(0,1) = mdl_center[1]; 
      eigen_avg_mdl(0,2) = mdl_center[2]; 
      eigen_avg_ref(0,0) = ref_center[0]; 
      eigen_avg_ref(0,1) = ref_center[1]; 
      eigen_avg_ref(0,2) = ref_center[2]; 
      Eigen::Matrix<double,1,3> translation = eigen_rot *
                                              (-eigen_avg_mdl.transpose()) + 
                                              eigen_avg_ref.transpose();
      transform(0, 3) = translation(0, 0);
      transform(1, 3) = translation(0, 1);
      transform(2, 3) = translation(0, 2);
      n_superposed = 2;
    } else {
      // the two cannot be superposed within specified distance threshold
      // => just set n_superposed to 1 and generate the same transformation
      // as in n_pos == 1
      transform = geom::Mat4::Identity();
      transform.PasteTranslation(ref_pos[0] - mdl_pos[0]);
      n_superposed = 1;
    }
    return;
  }

  if(window_size > n_pos) {
    window_size = n_pos;
  }

  Eigen::Matrix<double, Eigen::Dynamic, 3> eigen_mdl_pos = \
  Eigen::Matrix<double, Eigen::Dynamic, 3>::Zero(n_pos, 3);
  Eigen::Matrix<double, Eigen::Dynamic, 3> eigen_ref_pos = \
  Eigen::Matrix<double, Eigen::Dynamic, 3>::Zero(n_pos, 3);
  for(int i = 0; i < n_pos; ++i) {
  	eigen_mdl_pos(i, 0) = mdl_pos[i][0];
  	eigen_mdl_pos(i, 1) = mdl_pos[i][1];
  	eigen_mdl_pos(i, 2) = mdl_pos[i][2];
  	eigen_ref_pos(i, 0) = ref_pos[i][0];
  	eigen_ref_pos(i, 1) = ref_pos[i][1];
  	eigen_ref_pos(i, 2) = ref_pos[i][2];
  }

  std::vector<int> start_indices;
  int last_window_idx = n_pos - window_size;
  int n_windows = last_window_idx + 1;

  if(n_windows <= max_windows) {
    start_indices.resize(n_windows);
    for(int i = 0; i < n_windows; ++i) {
      start_indices[i] = i;
    }
  } else {
    start_indices.resize(max_windows);
    Real tmp = 0.0;
    Real delta = static_cast<Real>(last_window_idx) / (max_windows - 1);
    for(int i = 0; i < max_windows; ++i) {
      start_indices[i] = std::round(tmp);
      tmp += delta;
    }
  }
  
  size_t max_n = 0;
  Eigen::Matrix<double,1,3> eigen_avg_mdl = Eigen::Matrix<double,1,3>::Zero();
  Eigen::Matrix<double,1,3> eigen_avg_ref = Eigen::Matrix<double,1,3>::Zero();
  Eigen::Matrix<double,3,3> eigen_rotation = Eigen::Matrix<double,3,3>::Identity();
  Eigen::Matrix<double,1,3> eigen_avg_mdl_tmp = Eigen::Matrix<double,1,3>::Zero();
  Eigen::Matrix<double,1,3> eigen_avg_ref_tmp = Eigen::Matrix<double,1,3>::Zero();
  Eigen::Matrix<double,3,3> eigen_rotation_tmp = Eigen::Matrix<double,3,3>::Identity();
  
  for(int start_idx: start_indices) {

    std::vector<int> indices(window_size);
    for(int i = 0; i < window_size; ++i) {
      indices[i] = start_idx + i;
    }

    SuperposeIterative(eigen_mdl_pos, eigen_ref_pos,
                       10, distance_thresh, indices,
                       eigen_avg_mdl_tmp,
                       eigen_avg_ref_tmp,
                       eigen_rotation_tmp);

    if(indices.size() > max_n) {
      max_n = indices.size();
      eigen_avg_mdl = eigen_avg_mdl_tmp;
      eigen_avg_ref = eigen_avg_ref_tmp;
      eigen_rotation = eigen_rotation_tmp;
    }
  }

  // construct transform

  // there are three transformation to be applied to reach ref_pos from mdl_pos:
  // 1: shift mdl_pos to center
  // 2: apply estimated rotation
  // 3: shift onto average of ref_pos
  Eigen::Matrix<double,1,3> translation = eigen_rotation *
                                          (-eigen_avg_mdl.transpose()) + 
                                          eigen_avg_ref.transpose();


  transform = geom::Mat4();
  transform(0, 3) = translation(0, 0);
  transform(1, 3) = translation(0, 1);
  transform(2, 3) = translation(0, 2);
  transform(0, 0) = eigen_rotation(0, 0);
  transform(0, 1) = eigen_rotation(0, 1);
  transform(0, 2) = eigen_rotation(0, 2);
  transform(1, 0) = eigen_rotation(1, 0);
  transform(1, 1) = eigen_rotation(1, 1);
  transform(1, 2) = eigen_rotation(1, 2);
  transform(2, 0) = eigen_rotation(2, 0);
  transform(2, 1) = eigen_rotation(2, 1);
  transform(2, 2) = eigen_rotation(2, 2);

  n_superposed = max_n;
}

}}} // ns
