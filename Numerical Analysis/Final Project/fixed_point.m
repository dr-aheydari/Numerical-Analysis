 fp = @(wFixPoint) pi/2 - cos(wFixPoint);

                wFixPoint = 0;

                while (abs(wFixPoint - fp(wFixPoint))> 1e-9)

                    wFixPoint = fp(wFixPoint);


                end

                root = wFixPoint;
