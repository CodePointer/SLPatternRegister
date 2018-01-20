center_point = [0, 0, 35];
rad = 5.0;

cam_mat = [2426.231977104875, 0, 634.7636169096244;
 0, 2423.019176341463, 422.664603232326;
 0, 0, 1];
fx = cam_mat(1,1);
fy = cam_mat(2,2);
dx = cam_mat(1,3);
dy = cam_mat(2,3);
x_c = center_point(1);
y_c = center_point(2);
z_c = center_point(3);

for frm_idx = 0:20
  x_c = center_point(1) - frm_idx * 0.1;
  y_c = center_point(2);
  z_c = center_point(3);
  depth_mat = zeros(1024, 1280);
  dyna_mat = zeros(1024, 1280);
  for h = 1:1024
    for w = 1:1280
      M = (w - dx - 1) / fx;
      N = (h - dy - 1) / fy;
      A = M^2 + N^2 + 1;
      B = -2*M*x_c - 2*N*y_c - 2*z_c;
      C = x_c^2 + y_c^2 + z_c^2 - rad^2;
      if B^2 - 4*A*C < 0
        depth_mat(h, w) = -1;
        dyna_mat(h, w) = 15.0;
      else
        depth_mat(h, w) = (-B - sqrt(B^2-4*A*C)) / (2*A);
      end
    end
  end
end