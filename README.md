# read_pcl
     
## FPFH

https://pcl.readthedocs.io/projects/tutorials/en/latest/fpfh_estimation.html

worflow

for each point p in cloud P

  1. pass 1:

     1. get the nearest neighbors of :math:`p`

     2. for each pair of :math:`p, p_i` (where :math:`p_i` is a neighbor of :math:`p`, compute the three angular values

     3. bin all the results in an output SPFH histogram

  2. pass 2:

     1. get the nearest neighbors of :math:`p`

     3. use each SPFH of :math:`p` with a weighting scheme to assemble the FPFH of :math:`p`:
     
### FPFH sample
https://github.com/PointCloudLibrary/pcl/blob/master/features/include/pcl/features/fpfh.h
https://github.com/PointCloudLibrary/pcl/blob/master/examples/features/example_fast_point_feature_histograms.cpp

#### 继承关系
FPFHEstimation->FeatureFromNormals->Feature, Feature里面有公有方法compute计算Feature，内部调用的是虚函数computeFeature, FPFHEstimation实现了 虚函数computeFeature 

  pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_estimation;
  // Provide the original point cloud (without normals)
  
  fpfh_estimation.setInputCloud (cloud);
  
  // Provide the point cloud with normals

  fpfh_estimation.setInputNormals (cloud_with_normals);

  // fpfhEstimation.setInputWithNormals(cloud, cloudWithNormals); PFHEstimation does not have this function

  // Use the same KdTree from the normal estimation
  
  fpfh_estimation.setSearchMethod (tree);

  pcl::PointCloud<pcl::FPFHSignature33>::Ptr pfh_features (new pcl::PointCloud<pcl::FPFHSignature33>);

  fpfh_estimation.setRadiusSearch (0.2);

  // Actually compute the spin images

  fpfh_estimation.compute (*pfh_features);
  
### FPFHEstimation 

 
https://github.com/PointCloudLibrary/pcl/blob/master/features/include/pcl/features/fpfh.h

 
  template <typename PointInT, typename PointNT, typename PointOutT = pcl::FPFHSignature33>
  class FPFHEstimation : public FeatureFromNormals<PointInT, PointNT, PointOutT>  

feature 
https://github.com/PointCloudLibrary/pcl/blob/master/features/include/pcl/features/feature.h
  template <typename PointInT, typename PointOutT>
  class Feature : public PCLBase<PointInT>
  {
  pubic:
  ...
        void
      compute (PointCloudOut &output);//compute 调用虚函数computeFeature
      
   private:
       computeFeature (PointCloudOut &output) = 0;
  }
  
   template <typename PointInT, typename PointNT, typename PointOutT>
  class FeatureFromNormals : public Feature<PointInT, PointOutT>
  
  
https://github.com/PointCloudLibrary/pcl/blob/master/features/include/pcl/features/fpfh.h
https://github.com/PointCloudLibrary/pcl/blob/master/features/include/pcl/features/impl/feature.hpp

  template <typename PointInT, typename PointNT, typename PointOutT = pcl::FPFHSignature33>
  class FPFHEstimation : public FeatureFromNormals<PointInT, PointNT, PointOutT>

template <typename PointInT, typename PointOutT> void
Feature<PointInT, PointOutT>::compute (PointCloudOut &output)
{
  if (!initCompute ())
  {
    output.width = output.height = 0;
    output.clear ();
    return;
  }

  // Copy the header
  output.header = input_->header;

  // Resize the output dataset
  if (output.size () != indices_->size ())
    output.resize (indices_->size ());

  // Check if the output will be computed for all points or only a subset
  // If the input width or height are not set, set output width as size
  if (indices_->size () != input_->points.size () || input_->width * input_->height == 0)
  {
    output.width = indices_->size ();
    output.height = 1;
  }
  else
  {
    output.width = input_->width;
    output.height = input_->height;
  }
  output.is_dense = input_->is_dense;

  // Perform the actual feature computation
  computeFeature (output);//调用FPFH 实现

  deinitCompute ();
}

FPFH 实现
https://github.com/PointCloudLibrary/pcl/blob/master/features/include/pcl/features/impl/fpfh.hpp
template <typename PointInT, typename PointNT, typename PointOutT> void
pcl::FPFHEstimation<PointInT, PointNT, PointOutT>::computeFeature (PointCloudOut &output)


template <typename PointInT, typename PointNT, typename PointOutT> 
void pcl::FPFHEstimation<PointInT, PointNT, PointOutT>::computeFeature (PointCloudOut &output)
{
  // Allocate enough space to hold the NN search results
  // \note This resize is irrelevant for a radiusSearch ().
  pcl::Indices nn_indices (k_);
  std::vector<float> nn_dists (k_);

  std::vector<int> spfh_hist_lookup;
  computeSPFHSignatures (spfh_hist_lookup, hist_f1_, hist_f2_, hist_f3_);

  output.is_dense = true;
  // Save a few cycles by not checking every point for NaN/Inf values if the cloud is set to dense
  if (input_->is_dense)
  {
    // Iterate over the entire index vector
    for (std::size_t idx = 0; idx < indices_->size (); ++idx)
    {
      if (this->searchForNeighbors ((*indices_)[idx], search_parameter_, nn_indices, nn_dists) == 0)
      {
        for (Eigen::Index d = 0; d < fpfh_histogram_.size (); ++d)
          output[idx].histogram[d] = std::numeric_limits<float>::quiet_NaN ();

        output.is_dense = false;
        continue;
      }

      // ... and remap the nn_indices values so that they represent row indices in the spfh_hist_* matrices
      // instead of indices into surface_->points
      for (auto &nn_index : nn_indices)
        nn_index = spfh_hist_lookup[nn_index];

      // Compute the FPFH signature (i.e. compute a weighted combination of local SPFH signatures) ...
      weightPointSPFHSignature (hist_f1_, hist_f2_, hist_f3_, nn_indices, nn_dists, fpfh_histogram_);

      // ...and copy it into the output cloud
      std::copy_n(fpfh_histogram_.data (), fpfh_histogram_.size (), output[idx].histogram);
    }
  }
  else
  {
    // Iterate over the entire index vector
    for (std::size_t idx = 0; idx < indices_->size (); ++idx)
    {
      if (!isFinite ((*input_)[(*indices_)[idx]]) ||
          this->searchForNeighbors ((*indices_)[idx], search_parameter_, nn_indices, nn_dists) == 0)
      {
        for (Eigen::Index d = 0; d < fpfh_histogram_.size (); ++d)
          output[idx].histogram[d] = std::numeric_limits<float>::quiet_NaN ();

        output.is_dense = false;
        continue;
      }

      // ... and remap the nn_indices values so that they represent row indices in the spfh_hist_* matrices
      // instead of indices into surface_->points
      for (auto &nn_index : nn_indices)
        nn_index = spfh_hist_lookup[nn_index];

      // Compute the FPFH signature (i.e. compute a weighted combination of local SPFH signatures) ...
      weightPointSPFHSignature (hist_f1_, hist_f2_, hist_f3_, nn_indices, nn_dists, fpfh_histogram_);

      // ...and copy it into the output cloud
      std::copy_n(fpfh_histogram_.data (), fpfh_histogram_.size (), output[idx].histogram);
    }
  }
}

