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

% Set pattern
pat_mat = zeros(800, 1280);
for h = 1:3
  for w = 1:3
    coord = pro_point_coord{h, w};
    pat_mat(coord(1), coord(2)) = 1.0;
  end
end
imwrite(pat_mat, 'pattern_9p0.png');