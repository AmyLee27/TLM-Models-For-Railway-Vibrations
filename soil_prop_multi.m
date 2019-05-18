function [A_soil,B_soil,element_number,elements,element_coor,cinf] = soil_prop_multi(layer_heights, damping,var1_d,var2_d,var3_d,var1_v,var2_v,var3_v)

k1max = 10;
if(length(layer_heights) == 1 && layer_heights == 0) % If homogenous half-space
    var1_v = [var1_v var1_v]';               % replicate values to account for seperate supporting half-space
    var2_v = [var2_v var2_v]';
    var3_v = [var3_v var3_v]';
    damping = [damping damping]';
    layer_heights = [3 0]';                  % Set default value of upper TLM zone to 5m
    cinf = 1;                               % halfspace(1=yes,0=no)
elseif(length(layer_heights) == 1 && layer_heights ~= 0) % If homogenous rigid base
    cinf = 0;                               % halfspace(1=yes,0=no)
elseif(length(layer_heights) > 1 && layer_heights(end)==0) % If multi-layer with homogenous half-space
    cinf = 1;   % halfspace(1=yes,0=no)
%     layer_heights(end) = layer_heights(end-1)+8; 
   %layer_heights(end) = layer_heights(end-1)+1;
elseif(length(layer_heights) > 1 && layer_heights(end)~=0) % If multi-layer with rigid base
    cinf = 0;                               % halfspace(1=yes,0=no)
end

% Generate input var description matricies
total_soil_layers = length(damping);
var1_d = repmat(var1_d, [total_soil_layers 1]); % Vertical stack of row vectors
var2_d = repmat(var2_d, [total_soil_layers 1]);
var3_d = repmat(var3_d, [total_soil_layers 1]);

% Calculate elastodynamic properties for all layers
for x=1:length(damping)
    [E(x) G(x) K(x) lambda(x) M(x) poisson(x) density(x) Vs(x) Vp(x) Vr(x)]= elastodynamic(var1_d(x),var2_d(x),var3_d(x),var1_v(x),var2_v(x),var3_v(x));    
end

% Define stability criteria
 max_TL_per_wavelen = 8;
% calculate the wavelength
wavelength = 2*pi/k1max; % wavelength(lamda)= speed/frequency (c/f) = 2*pi/wavenumber(2*pi/k)
% determine the thickness of thin layer and the total number of thin layers used in the model.
thinlayer_thickness = wavelength/max_TL_per_wavelen;  % thin layer thickness = wavelength/8 (usually 8 pr 10)

% calculate the number of thin layers in each layer and label them
TP = 0;
hh = [0 layer_heights'];
element_coor = zeros(1000000,1);  % define a large matrix for it first
element_coor(1,1) = 0;  % set up the first element to 0

% Define nodal coordinates etc whether a half-space or not
if(layer_heights(end)==0)               % Half-space
    for i = 1:(total_soil_layers-1)
        kh(i) = hh(i+1) - hh(i);
        p(i) = ceil(kh(i)/thinlayer_thickness);
        TK = TP + 1;
        TP = p(i) + TP;
        adjust_thickness(i) = kh(i)/p(i); % calculate the adjusted thin-layer thickness
        
        for j = TK:TP
            element_number(j,:) = [j 2*j-1 2*j 2*j+1]; % Matrix of elements/nodes number
            A_soil(j,:) = [j E(i) poisson(i) damping(i) density(i) j];
            
            for k = (2*j):(2*j+1)
                element_coor(k,1) = element_coor(k-1,1) + adjust_thickness(i)/2;
                % calculate the nodes positions in the thin-layers
            end
        end
        
    end
    B_soil = [E(end),poisson(end),damping(end),density(end)];
    element_coor((2*TP+2):end) = []; % Remove all the zero rows to truncate the matrix into a proper size
    elements = (1:(2*TP+1))';  % list the nodes number
else            % Rigid base
    for i = 1:(total_soil_layers)
        kh(i) = hh(i+1) - hh(i);
        p(i) = ceil(kh(i)/thinlayer_thickness);
        TK = TP + 1;
        TP = p(i) + TP;
        adjust_thickness(i) = kh(i)/p(i); % calculate the adjusted thin-layer thickness
        
        for j = TK:TP
            element_number(j,:) = [j 2*j-1 2*j 2*j+1]; % Matrix of elements/nodes number
            A_soil(j,:) = [j E(i) poisson(i) damping(i) density(i) j];
            
            for k = (2*j):(2*j+1)
                element_coor(k,1) = element_coor(k-1,1) + adjust_thickness(i)/2;
                % calculate the nodes positions in the thin-layers
            end
        end
        
    end
    B_soil = [E(end),poisson(end),damping(end),density(end)];
    element_coor((2*TP+2):end) = []; % Remove all the zero rows to truncate the matrix into a proper size
    elements = (1:(2*TP+1))';  % list the nodes number
end
