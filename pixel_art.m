% Project v1

% Initialization
% Input image
file_name = 'trees.tif';
[aind,amap]	=   imread(file_name ,'TIF');
image_input = ind2rgb(aind,amap);
image_input = rgb2lab(image_input);
[h_in, w_in, nb_color] = size(image_input);
M = w_in * h_in;

% Output image
w_out = w_in / 2;
h_out = h_in / 2;
N = w_out*h_out;
image_output = zeros(h_in,w_in,3); 
% [L,numLabels] = superpixel(image_input,N); % not he exact method, as the
% one explained from the paper is slightly different. Cannot process double

super = zeros(h_out, w_out,3);
for k = 1:h_out
    for l = 1:w_out
        super(k,l,1) = k*(10+1)/2; % average position : row
        super(k,l,2) = l*(10+1)/2; % average position : column
    end
end


% Colors
max_color = 16; % max number of color of the output image 
mean_color = zeros(1,3); % mean color of the input image
for k = 1:nb_color
    mean_color(1) = mean(image_input(:,:,1),[1,2]);
    mean_color(2) = mean(image_input(:,:,2),[1,2]);
    mean_color(3) = mean(image_input(:,:,3),[1,2]);
end
palette = zeros(max_color,3);
palette(1,1) = mean_color(1);
palette(1,2) = mean_color(2);
palette(1,3) = mean_color(3);

% all superpixels starts with the first color of the palette
for k = 1:h_out
    for l = 1:w_out
        super(k,l,3) = 1; % position of the color in palette
    end
end


% Critical temperature
if w_out > h_out
    Tc = 2*var(image_input,0,1);
else
    Tc = 2*var(image_input,0,2);
end

T = 1.1*Tc;

m = 45; % varaible depending on the case

Tf = T-1; % temporary, for 1 loop

% While loop
while T >Tf
    % for all pixel of the input image
    for k = 1:h_in
        for l = 1:w_in
            d = zeros(h_out,w_out);
            % compare with all superpixel
            for a = 1:h_out
                for b = 1:w_out
                    difference_color_bewteen_pixel_and_super_pixel = abs(image_input(k, l, 1) - palette(super(a,b,3),1)) + abs(image_input(k, l, 2) - palette(super(a,b,3),2)) + abs(image_input(k, l, 3) - palette(super(a,b,3),3));
                    difference_position_bewteen_pixel_and_super_pixel = abs(k-super(a,b,1)) + abs(l-super(a,b,2));
                    d(a,b) = difference_color_bewteen_pixel_and_super_pixel + m * sqrt(N/M)* difference_position_bewteen_pixel_and_super_pixel;
                end
            end           
            % take the min of d and associate the pixel with the superpixel
            [row,col] = find(d==min(d(:)));
            image_output(k,l,1) = super(row,col,1);
            image_output(k,l,2) = super(row,col,2);
            image_output(k,l,3) = super(row,col,3);
        end
    end
   
T= Tf; % temporary, for 1 loop

end % end of the while loop