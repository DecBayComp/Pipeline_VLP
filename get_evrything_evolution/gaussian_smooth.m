function vecteur = gaussian_smooth(vecteur, n)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/
windowWidth              = int32(n);
halfWidth                = windowWidth / 2;
gaussFilter              = gausswin(n);
gaussFilter              = gaussFilter / sum(gaussFilter); %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%/

n_vecteur      = length(vecteur);
vecteur        = [vecteur(1)*ones(halfWidth,1);vecteur;vecteur(n_vecteur)*ones(halfWidth,1)];
vecteur        = conv(vecteur, gaussFilter, 'same');
vecteur        = vecteur(halfWidth+1:length(vecteur)-halfWidth);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
