function save_intensity_16bit(intensity, filename)

scaled_intensity = intensity ./ max(intensity(:));
scaled_intensity = uint16((2^16 - 1) * scaled_intensity);
imwrite(scaled_intensity(:, :, 1), filename)
for i=2:size(intensity, 3)
    imwrite(scaled_intensity(:, :, i), filename, 'WriteMode', 'append')
end

end
