function [LGFPang, LGFPmag]=localFrissRadPatt(Ig)
dm = size(Ig,3);
if dm >1
    Ig = double(rgb2gray(Ig));
else
    Ig = double(Ig);
end

if ~isa(Ig, 'double')
    Ig = im2double(Ig);
end

R = 1;
P = 8;
[rw, cl] = size(Ig);
LGFPmag = zeros(rw,cl);
LGFPang = zeros(rw,cl);
k = R;
ang = 2*pi/P;
I = Ig;
for i = 1+k:rw-k
    for j = 1+k:cl-k
        ForceX = 0;
        ForceY = 0;
        for p = 0:P-1
            x = -(R(1)*sin((p)*ang));
            if x>0
                x = x+.1;
            else
                x = x-.1;
            end
            y = (R(1)*cos((p)*ang));
            if y>0
                y = y+.1;
            else
                y = y-.1;
            end
            r = round(i+x);
            c = round(j+y);
            lamb = 650-(I(i,j)*0.72);
            ForceX = ForceX+((cos(rem(p*ang+pi,2*pi))*I(r,c)*(4*pi*((r-i)^2+(c-j)^2))^2)/(1*I(i,j)*1^2));
            ForceY = ForceY+((sin(rem(p*ang+pi,2*pi))*I(r,c)*(4*pi*((r-i)^2+(c-j)^2))^2)/(1*I(i,j)*1^2));
        end
        LGFPmag(i,j) = sqrt(ForceX^2+ForceY^2);
        LGFPang(i,j) = (atan2(ForceY,ForceX)+pi)*180/pi;
    end
end