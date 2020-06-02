function [rotatedallaxes, integerdirections] = rotateaxes(allaxes,R)
%rotateaxes - rotate axes%
%Inputs%
%allaxes - n by 4 by 3 matrix of axis vectors and lengths from Fit_principal axes function%
%R - rotation matrix from transform axes function, rotates tip system
%to crystal system% crystal = R*tip (column vectors).%
%Outputs%
%rotatedallaxes, same matrix as all axes but where each vector has been
%rotated, all still with corresponding lengths
%intergerdirections, attempt to out put classic crystal inter directions
%closest to vectors%
rotatedallaxes = allaxes;
for i = 1:size(allaxes,1)
    extract(1:3,1:3) = allaxes(i,1:3,:); %extacrt matrix of column vectors for each axis
    rotatedallaxes(i,1:3,:) = R*extract; %rotate vectors to the crystal basis%
    
end
integerdirections = rotatedallaxes;
for i = 1:size(allaxes)
    for j = 1:3
        x = integerdirections(i,1,j);
        y = integerdirections(i,2,j);
        z = integerdirections(i,3,j);
        p = 1000; %number of factors being tried%
        t = 0.15; %tolerance%
        factor = zeros(p,1);
        for k = 1:p 
            f = k/p*1000; %current maximum factor tried is 100%
            if (f*x - round(f*x,0) > -t) && (f*x - round(f*x,0) < t)...
                    && (f*y - round(f*y,0) > -t) && (f*y - round(f*y,0) < t)...
                    && (f*z - round(f*z,0) > -t) && (f*z - round(f*z,0) < t)
                factor(k,1) = f;
            else
                factor(k,1) = NaN;
            end
        end
        A = integerdirections(i,1:3,j);
        integerdirections(i,1:3,j) = round(A*min(factor));
    end
end
end

