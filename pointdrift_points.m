function [cents] = pointdrift_points(im) 
    config = get_config();
    orig = im;
    
    if config('pointdrift_use_dog')
        im = pointdrift_dog(im);
    else
        % FIXME constant size here
        %size_threshold = 40;

        % smooth the image
        %im = histeq(im);
        im = imtophat(im, strel('disk', config('pointdrift_tophat_radius')));

        if config('pointdrift_gaussian')
            h = fspecial('gaussian');
            im = imfilter(im, h);
        end
    end

    if config('pointdrift_quantize')
        % quantize into 3 classes
        numclasses = config('pointdrift_quantize_classes');
        quant = imquantize(im, multithresh(im, numclasses-1));
        % grab class 3, which is (presumably) the brightest points
        bw = quant==numclasses;
        if config('debug') ; figure,imshow(label2rgb(quant, 'lines', 'k'));title('quant'); end
    else
        bw = im2bw(im);
    end
    
    %bw = imdilate(bw, strel('disk',1));
    %rm = bwareaopen(bw, size_threshold);
    %bw(rm) = 0;
    
    %bw = im2bw(b, graythresh(im));
    %figure,imshow(bw);
    if config('debug') ; figure,imshow(bw),title('bw'); end

    
    props = regionprops(bw, 'Centroid');
    cents = cat(1, props.Centroid);
    cents = cents(isfinite(cents(:, 1)), :);
    %figure,imshow(im),hold on,cpd_plot_iter(cents, cents);
    if config('debug')
        figure,imshow(orig),hold on,scatter(cents(:, 1), cents(:, 2), 15, repmat([0,0,1] ,size(cents,1),1), '+', 'blue');title('points');
    end
end

function [dogim] = pointdrift_dog(im) 
    config = get_config();
    dog = fspecial('Gaussian', 21, 15) - fspecial('Gaussian', 21, 20);
    dogim = conv2(double(im), dog, 'same');
    if config('debug') ; figure,imshow(dogim);title('dog'); end
end