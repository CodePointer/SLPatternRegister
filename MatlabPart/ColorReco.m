%% Prepare
point_range = [300, 500; 700, 1000];
% Generate 9 points from point_range;
pro_point_coord = cell(3, 3);
pro_point_coord{1, 1} = [point_range(1, 1); point_range(1, 2)];
pro_point_coord{1, 3} = [point_range(1, 1); point_range(2, 2)];
pro_point_coord{3, 1} = [point_range(2, 1); point_range(1, 2)];
pro_point_coord{3, 3} = [point_range(2, 1); point_range(2, 2)];
pro_point_coord{1, 2} = 0.5 * (pro_point_coord{1,1} + pro_point_coord{1,3});
pro_point_coord{2, 1} = 0.5 * (pro_point_coord{1,1} + pro_point_coord{3,1});
pro_point_coord{3, 2} = 0.5 * (pro_point_coord{3,1} + pro_point_coord{3,3});
pro_point_coord{2, 3} = 0.5 * (pro_point_coord{1,3} + pro_point_coord{3,3});
pro_point_coord{2, 2} = 0.5 * (pro_point_coord{1,2} + pro_point_coord{3,2});
% Match filter
match_filter = fspecial('gauss', 5, 1.5);

%% Search For 9 points from point_range
% Read
img_9p = imread('img_9p.png');
img_9p = rgb2gray(img_9p);
[kImgHeight, kImgWidth] = size(img_9p);
% img_filtered = imfilter(img_9p, match_filter);
img_filtered = img_9p;

% Find max points
max_points = [];
for h = 5:kImgHeight-5
  for w = 5:kImgWidth-5
    if (img_filtered(h, w) == max(max(img_filtered(h-2:h+2, w-2:w+2)))) ...
        && (img_filtered(h, w) ~= min(min(img_filtered(h-2:h+2, w-2:w+2))))
      max_points = [max_points, [h; w]];
    end
  end
end

% From mahhaton find coords
cam_point_coord = cell(3, 3);
mah_dist = zeros(9, 1);
for i = 1:9
  mah_dist(i) = max_points(1, i) + max_points(2, i);
end
[~, sorted_idx] = sort(mah_dist);
cam_point_coord{1, 1} = max_points(:, sorted_idx(1));
sub_set_points = max_points(:, sorted_idx(2:3));
[~, sub_idx] = sort(sub_set_points(1,:));
cam_point_coord{1, 2} = sub_set_points(:, sub_idx(1));
cam_point_coord{2, 1} = sub_set_points(:, sub_idx(2));
sub_set_points = max_points(:, sorted_idx(4:6));
[~, sub_idx] = sort(sub_set_points(1,:));
cam_point_coord{1, 3} = sub_set_points(:, sub_idx(1));
cam_point_coord{2, 2} = sub_set_points(:, sub_idx(2));
cam_point_coord{3, 1} = sub_set_points(:, sub_idx(3));
sub_set_points = max_points(:, sorted_idx(7:8));
[~, sub_idx] = sort(sub_set_points(1,:));
cam_point_coord{2, 3} = sub_set_points(:, sub_idx(1));
cam_point_coord{3, 2} = sub_set_points(:, sub_idx(2));
cam_point_coord{3, 3} = max_points(:, sorted_idx(9));

%% Calculate transfer mat
pro_points = zeros(9, 2);
cam_points = zeros(9, 2);
for h = 1:3
  for w = 1:3
    i = (h-1)*3+w;
    pro_points(i, :) = pro_point_coord{h, w};
    cam_points(i, :) = cam_point_coord{h, w};
  end
end
tform = fitgeotrans(cam_points, pro_points, 'projective');
% for h_i = 1:16
%     for w_i = 1:16
%         tmp_pos = [h_i, w_i, 1] * tform.T;
%         h_image = round(tmp_pos(1) / tmp_pos(3));
%         w_image = round(tmp_pos(2) / tmp_pos(3));
%         position_set{h_i, w_i} = [h_image; w_image];
%     end
% end

%% Set x_pro, y_pro mat
x_pro_mat = - ones(kImgHeight, kImgWidth);
y_pro_mat = - ones(kImgHeight, kImgWidth);
% mask_mat = zeros(kImgHeight, kImgWidth);
img_obs = imread('img_obs.png');
mask_mat = img_obs > 20;
for h = 1:kImgHeight
  for w = 1:kImgWidth
    if mask_mat(h, w) 
      tmp_pos = [h, w, 1] * tform.T;
      y_pro_mat(h, w) = tmp_pos(2) / tmp_pos(3) - 1;
      x_pro_mat(h, w) = tmp_pos(1) / tmp_pos(3) - 1;
    end
  end
end