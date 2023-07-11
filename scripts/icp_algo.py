
import numpy as np
import
class mag_icp:
    """Class for the MagNav ICP algorithm."""

    def __init__(self):
        """Attributes and initializers for the mag_icp class."""
        self.partial_update





    def icp_algo(self, model, data, true_path, initial_path, est_path, maxIter=40, minIter=5, critFun=0, thres=1e-5):
        """
        Python implementation of the ICP algorithm from MATLAB.
        RETURNS: TR, TT, data, res
        Reference:
        Bergström, P. and Edlund, O. 2014, 'Robust registration of point sets using iteratively reweighted least squares'
        Computational Optimization and Applications, vol 58, no. 3, pp. 543-561, 10.1007/s10589-014-9643-2
        """

        if not model or if not data:
            print("Something is wrong with the model and/or data points.")


        # TODO: Convert to Python Code
        #Size of model points and data points
        # if (size(model, 2) < size(model, 1))
        #     mTranspose = true;
        #     m = size(model, 2);
        #     M = size(model, 1);
        # else
        #     mTranspose = false;
        #     m = size(model, 1);
        #     M = size(model, 2);
        # end
        #
        # if (size(data, 2) < size(data, 1))
        #     data = data';
        # end
        #
        # if m~=size(data, 1)
        # error('The dimension of the model points and data points must be equal');
        # end
        #
        # N = size(data, 2);



        # TODO: Convert to Python Code
        # Create closest point search structure
        # if m < 4
        #     if mTranspose
        #         disp("Delaudnay Triangulation, mTranspose True")
        #         DT = delaunayTriangulation(model);
        #     else
        #         disp("Delaudnay Triangulation, mTranspose False")
        #         DT = delaunayTriangulation(model');
        #     end
        # else
        #     DT = [];
        #     resid = zeros(N, 1);
        #     vi = ones(N, 1);
        # end



        # TODO: Convert to Python Code
        # # Initiate weights (Only for robust criterion)
        # if critFun > 0
        #     wghs = ones(N, 1);
        # end
        #
        # # Initiate transformation
        # TR = eye(m);
        # TT = zeros(m, 1);



        # TODO: Convert to Python Code
        # Start the ICP algorithm
        # res = 9e99;
        #
        # for iter=1:maxIter
        #
        # oldres = res;
        #
        # %Find
        # closest
        # model
        # points
        # to
        # data
        # points
        # if isempty(DT)
        #     if mTranspose
        #         for i=1:N
        #             mival = 9e99;
        #             for j=1:M
        #                 val = norm(data(:, i)-model(j,:)');
        #                 if val < mival
        #                     mival = val;
        #                     vi(i) = j;
        #                     resid(i) = val;
        #                 end
        #             end
        #         end
        #     else
        #         for i=1:N
        #             mival = 9e99;
        #             for j=1:M
        #                 val = norm(data(:, i)-model(:, j));
        #                 if val < mival
        #                     mival=val;
        #                     vi(i)=j;
        #                     resid(i)=val;
        #                 end
        #             end
        #         end
        #     end
        # else
        #     [vi, resid] = nearestNeighbor(DT, data');
        # end







