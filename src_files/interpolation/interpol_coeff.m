function[ap,an,b] = interpol_coeff(faces,centroids,neighb,neighb_pos,phi,grad_phi,bd,S,neighb_bound,phi_bound,F,Schemes)

grad_phi_p = grad_phi;
grad_phi_n = grad_phi(neighb_pos,:);

grad_phi_f = grad_interp(faces,centroids,neighb,neighb_pos,phi,grad_phi,bd,S,neighb_bound,phi_bound);

pos = F>=0; % pos=1 implies upwind node is centroid

phi_neighb = phi(neighb_pos);

phi_c = phi;phi_c(~pos) = phi_neighb(~pos); % upwind of face
phi_d = phi_neighb;phi_d(~pos) = phi(~pos); % downwind of face

grad_phi_c = grad_phi_p;grad_phi_c(~pos,:) = grad_phi_n(~pos,:); % upwind of face
grad_phi_d = grad_phi_n;grad_phi_d(~pos,:) = grad_phi_p(~pos,:); % downwind of face

dcf = (faces - centroids);
dcf(~pos,:) = (faces(~pos,:) - neighb(~pos,:));

% Refer page 409 on uFVM book
% constants for different schemes

div_scheme = Schemes.div_scheme;

switch(div_scheme)
    case 'UD'
        c1 = 0; c2 = 0; % UD
    case 'CD'
        c1 = 0; c2 = 1; % CD
    case 'SOU'
        c1 = 2; c2 = -1; % SOU
    case 'FROMM'
        c1 = 1; c2 = 0; % FROMM
    case 'QUICK'
        c1 = 0.5; c2 = 0.5; % QUICK
end

phi_f = phi_c + dot(c1*grad_phi_c + c2*grad_phi_f,dcf,2);

ap = zeros(size(phi)); ap(pos) = 1;
an = zeros(size(phi)); an(~pos) =  1;
b  = phi_f - phi_c;


end