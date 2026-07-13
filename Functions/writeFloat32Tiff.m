function writeFloat32Tiff(A, path)
    A = single(A);
    t = Tiff(path, 'w');
    tag.ImageLength = size(A,1);
    tag.ImageWidth = size(A,2);
    tag.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tag.BitsPerSample = 32;
    tag.SamplesPerPixel = 1;
    tag.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tag.Photometric = Tiff.Photometric.MinIsBlack;
    for k = 1:size(A,3)
        setTag(t, tag); write(t, A(:,:,k));
        if k < size(A,3), writeDirectory(t); end
    end
    close(t);
end