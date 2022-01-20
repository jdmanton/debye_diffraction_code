function [pupil, aberration] = apply_aberration(pupil, defocus, ast_oblique, ast_vertical, coma_vertical, coma_horizontal, trefoil_vertical, trefoil_oblique, spherical_primary)

z_defocus = defocus .* (pupil.r.^2 - 1);
z_ast_oblique = ast_oblique .* pupil.r.^2 .* sind(2 .* pupil.phi);
z_ast_vertical = ast_vertical .* pupil.r.^2 .* cosd(2 .* pupil.phi);
z_coma_vertical = coma_vertical .* (3 .* pupil.r.^3 - 2 .* pupil.r) .* sind(pupil.phi);
z_coma_horizontal = coma_horizontal .* (3 .* pupil.r.^3 - 2 .* pupil.r) .* cosd(pupil.phi);
z_trefoil_vertical = trefoil_vertical .* pupil.r.^3 .* sind(3 .* pupil.phi);
z_trefoil_oblique = trefoil_oblique .* pupil.r.^3 .* cosd(3 .* pupil.phi);
z_spherical_primary = spherical_primary .* (6 .* pupil.r.^4 - 6 .* pupil.r.^2 + 1);

aberration = z_defocus + z_ast_oblique + z_ast_vertical + z_coma_vertical + z_coma_horizontal + z_trefoil_vertical + z_trefoil_oblique + z_spherical_primary;

pupil.phase = exp(1i .* aberration);

end
