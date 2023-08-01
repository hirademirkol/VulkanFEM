#include "MatrixFreeSparse.hpp"

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

namespace Eigen {
namespace internal {
  // template<typename scalar, int options, typename storageIndex, typename Rhs>
  // struct generic_product_impl<MatrixFreeSparse<scalar, options, storageIndex>, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  // : generic_product_impl_base<MatrixFreeSparse<scalar, options, storageIndex>,Rhs,generic_product_impl<MatrixFreeSparse<scalar, options, storageIndex>,Rhs> >
  // {
  //   typedef typename Product<MatrixFreeSparse<scalar, options, storageIndex>,Rhs>::Scalar Scalar;
 
  //   template<typename Dest>
  //   static void scaleAndAddTo(Dest& dst, const MatrixFreeSparse<scalar, options, storageIndex>& lhs, const Rhs& rhs, const Scalar& alpha)
  //   {
  //     // This method should implement "dst += alpha * lhs * rhs" inplace,
  //     // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
  //     assert(alpha==Scalar(1) && "scaling is not implemented");
  //     EIGEN_ONLY_USED_FOR_DEBUG(alpha);
 
  //     // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
  //     // but let's do something fancier (and less efficient):
  //     for(auto line : lhs.elementToNode)
  //     {
  //       for(auto&& [i , node1] : line )
  //         for(auto&& [j, node2] : line )
  //         {
  //           if(node1 == -1 || node2 == -1)
	// 				  continue;

  //           for(int c1 = 0; c1 < 3; c1++)
  //           {
  //             for(int c2 = 0; c2 < 3; c2++)
  //             {
  //               int iMatrix = node1*3 + c1;
  //               int jMatrix = node2*3 + c2;
  //               int elementIndex = i <= j ? (i*(47-i))/2 + j + 1 : (j*(47-j))/2 + i + 1;

  //               dst(i) += lhs.elementStiffnessMat(elementIndex) * rhs(j);
  //             }
  //           }
  //         }
  //     }
  //   };
  // };

  // // template struct generic_product_impl<MatrixFreeSparse<double, 0, int>, Eigen::Matrix<float,-1,1,0,-1,1>, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector

  template<typename Rhs>
  struct generic_product_impl<MatrixFreeSparse, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
  : generic_product_impl_base<MatrixFreeSparse,Rhs,generic_product_impl<MatrixFreeSparse,Rhs> >
  {
    typedef typename Product<MatrixFreeSparse,Rhs>::Scalar Scalar;
 
    template<typename Dest>
    static void scaleAndAddTo(Dest& dst, const MatrixFreeSparse& lhs, const Rhs& rhs, const Scalar& alpha)
    {
      // This method should implement "dst += alpha * lhs * rhs" inplace,
      // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
      assert(alpha==Scalar(1) && "scaling is not implemented");
      EIGEN_ONLY_USED_FOR_DEBUG(alpha);
 
      // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
      // but let's do something fancier (and less efficient):
      for(auto line : lhs.elementToNode)
      {
        for(auto&& [i , node1] : line )
          for(auto&& [j, node2] : line )
          {
            if(node1 == -1 || node2 == -1)
					  continue;

            for(int c1 = 0; c1 < 3; c1++)
            {
              for(int c2 = 0; c2 < 3; c2++)
              {
                int iMatrix = node1*3 + c1;
                int jMatrix = node2*3 + c2;
                int elementIndex = i <= j ? (i*(47-i))/2 + j + 1 : (j*(47-j))/2 + i + 1;

                dst(i) += lhs.elementStiffnessMat(elementIndex) * rhs(j);
              }
            }
          }
      }
    };
  };
  
}
}