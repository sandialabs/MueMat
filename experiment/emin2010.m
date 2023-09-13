% Generates tables of emin2010 paper

function emin2010

tic

ElasticityTest(2, 40:40:320, 1);
ElasticityTest(2, 40:40:320, 10);
ElasticityTest(2, 40:40:320, 100);

ElasticityTest(3, 5:5:40, 1);
ElasticityTest(3, 5:5:40, 10);
ElasticityTest(3, 5:5:40, 100);

toc

if 0
% extra
tic
ElasticityTest(3, 45:5:50, 1);
ElasticityTest(3, 45:5:50, 10);
ElasticityTest(3, 45:5:50, 100);
toc
end

ElasticityTable();

end
