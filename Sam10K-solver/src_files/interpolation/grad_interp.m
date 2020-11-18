function[grad_phi_f] = grad_interp(faces,centroids,neighb,neighb_pos,phi,grad_phi,bd,S,neighb_bound,phi_bound)

phi_n = phi(neighb_pos);
phi_p = phi;
grad_phi_x = grad_phi(:,1);
grad_phi_y = grad_phi(:,2);
fx = vecnorm(faces-neighb,2,2)./vecnorm(centroids-neighb,2,2);
fx(abs(fx)==Inf)=0;
fx(isnan(fx)) = 0;
grad_phi_f = [grad_phi_x.*fx + (1-fx).*grad_phi_x(neighb_pos),...
                grad_phi_y.*fx + (1-fx).*grad_phi_y(neighb_pos)];
bd = logical(bd);
d = (neighb-centroids);            
db = (faces - centroids);
dn = dot(db,S,2).*S./(S(:,1).^2 + S(:,2).^2);
phi_n(bd,:)=  phi_bound(neighb_bound(bd));
d(bd,:)=dn(bd,:);
e = d./vecnorm(d,2,2);

grad_phi_f(bd,:) = 0;
% Refer page 289 of uFVM book
grad_phi_f = [grad_phi_f(:,1) + ((phi_n-phi_p)./vecnorm(d,2,2) - dot(grad_phi_f,e,2)).*e(:,1),...
          grad_phi_f(:,2) + ((phi_n-phi_p)./vecnorm(d,2,2) - dot(grad_phi_f,e,2)).*e(:,2)];


end