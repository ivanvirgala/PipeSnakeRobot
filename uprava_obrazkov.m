% uprava obrazkov....z jpg obrazku odstrani nadbytocne biele miesta aby
% nebol obrazok roztiahnuty zbytocne
image = imread('Fig/traveledDistanceCOM_viscousVScoulomb2.jpg');
gray_image = rgb2gray(image);
non_white_columns = any(gray_image < 255, 1);
first_non_white_column = find(non_white_columns, 1, 'first');
last_non_white_column = find(non_white_columns, 1, 'last');
cropped_image = image(:, first_non_white_column:last_non_white_column, :);
imwrite(cropped_image, 'Fig/traveledDistanceCOM_viscousVScoulomb2Modified.jpg');