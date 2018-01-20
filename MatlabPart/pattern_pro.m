ori_pattern = imread('pattern_2size2color8P0.png');
d_pattern = double(ori_pattern) / 255.0;
set_pattern = d_pattern * 100 + 20;
h = fspecial('gauss', 7, 1.2);
final_pattern = imfilter(set_pattern, h);
imshow(final_pattern, [0, 255]);
save('pattern_gauss.txt', 'final_pattern', '-ascii')
imwrite(uint8(final_pattern), 'pattern_gauss.png');