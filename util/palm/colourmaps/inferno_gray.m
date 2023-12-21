function map = inferno_gray(n)
% Returns the "inferno" colourmap.

if nargin < 1
    n = size(get(gcf,'Colormap'),1);
end

values = [
   0.498039215686275   0.498039215686275   0.498039215686275
   0.490196078431373   0.490196078431373   0.490196078431373
   0.482352941176471   0.482352941176471   0.482352941176471
   0.474509803921569   0.474509803921569   0.474509803921569
   0.462745098039216   0.462745098039216   0.462745098039216
   0.454901960784314   0.454901960784314   0.454901960784314
   0.447058823529412   0.447058823529412   0.450980392156863
   0.439215686274510   0.439215686274510   0.443137254901961
   0.431372549019608   0.431372549019608   0.435294117647059
   0.423529411764706   0.423529411764706   0.427450980392157
   0.415686274509804   0.415686274509804   0.423529411764706
   0.407843137254902   0.407843137254902   0.415686274509804
   0.403921568627451   0.400000000000000   0.411764705882353
   0.396078431372549   0.392156862745098   0.403921568627451
   0.388235294117647   0.388235294117647   0.396078431372549
   0.380392156862745   0.380392156862745   0.392156862745098
   0.372549019607843   0.372549019607843   0.384313725490196
   0.364705882352941   0.364705882352941   0.380392156862745
   0.360784313725490   0.356862745098039   0.376470588235294
   0.352941176470588   0.349019607843137   0.368627450980392
   0.345098039215686   0.341176470588235   0.364705882352941
   0.341176470588235   0.341176470588235   0.364705882352941
   0.337254901960784   0.333333333333333   0.356862745098039
   0.329411764705882   0.325490196078431   0.352941176470588
   0.321568627450980   0.317647058823529   0.349019607843137
   0.317647058823529   0.309803921568627   0.345098039215686
   0.313725490196078   0.305882352941176   0.345098039215686
   0.309803921568627   0.301960784313725   0.341176470588235
   0.301960784313725   0.294117647058824   0.337254901960784
   0.301960784313725   0.290196078431373   0.337254901960784
   0.294117647058824   0.282352941176471   0.333333333333333
   0.290196078431373   0.274509803921569   0.329411764705882
   0.286274509803922   0.270588235294118   0.329411764705882
   0.282352941176471   0.266666666666667   0.325490196078431
   0.274509803921569   0.258823529411765   0.321568627450980
   0.274509803921569   0.254901960784314   0.321568627450980
   0.270588235294118   0.247058823529412   0.321568627450980
   0.270588235294118   0.243137254901961   0.321568627450980
   0.262745098039216   0.235294117647059   0.317647058823529
   0.262745098039216   0.231372549019608   0.321568627450980
   0.258823529411765   0.223529411764706   0.317647058823529
   0.258823529411765   0.219607843137255   0.317647058823529
   0.258823529411765   0.215686274509804   0.321568627450980
   0.254901960784314   0.207843137254902   0.317647058823529
   0.250980392156863   0.203921568627451   0.317647058823529
   0.247058823529412   0.196078431372549   0.317647058823529
   0.247058823529412   0.192156862745098   0.317647058823529
   0.250980392156863   0.188235294117647   0.321568627450980
   0.243137254901961   0.180392156862745   0.317647058823529
   0.247058823529412   0.176470588235294   0.317647058823529
   0.247058823529412   0.172549019607843   0.321568627450980
   0.247058823529412   0.172549019607843   0.321568627450980
   0.243137254901961   0.160784313725490   0.321568627450980
   0.243137254901961   0.156862745098039   0.321568627450980
   0.243137254901961   0.152941176470588   0.321568627450980
   0.247058823529412   0.149019607843137   0.325490196078431
   0.247058823529412   0.145098039215686   0.325490196078431
   0.247058823529412   0.145098039215686   0.325490196078431
   0.247058823529412   0.141176470588235   0.325490196078431
   0.247058823529412   0.133333333333333   0.325490196078431
   0.247058823529412   0.129411764705882   0.325490196078431
   0.250980392156863   0.125490196078431   0.325490196078431
   0.250980392156863   0.121568627450980   0.325490196078431
   0.250980392156863   0.121568627450980   0.329411764705882
   0.254901960784314   0.117647058823529   0.329411764705882
   0.254901960784314   0.113725490196078   0.329411764705882
   0.258823529411765   0.113725490196078   0.329411764705882
   0.258823529411765   0.109803921568627   0.329411764705882
   0.262745098039216   0.113725490196078   0.337254901960784
   0.266666666666667   0.109803921568627   0.337254901960784
   0.266666666666667   0.105882352941176   0.333333333333333
   0.270588235294118   0.105882352941176   0.337254901960784
   0.274509803921569   0.101960784313725   0.337254901960784
   0.274509803921569   0.101960784313725   0.337254901960784
   0.278431372549020   0.098039215686275   0.337254901960784
   0.286274509803922   0.101960784313725   0.341176470588235
   0.286274509803922   0.098039215686275   0.345098039215686
   0.290196078431373   0.098039215686275   0.345098039215686
   0.294117647058824   0.098039215686275   0.345098039215686
   0.298039215686275   0.094117647058824   0.345098039215686
   0.301960784313725   0.094117647058824   0.349019607843137
   0.309803921568627   0.098039215686275   0.352941176470588
   0.313725490196078   0.094117647058824   0.349019607843137
   0.313725490196078   0.094117647058824   0.349019607843137
   0.321568627450980   0.094117647058824   0.356862745098039
   0.329411764705882   0.094117647058824   0.352941176470588
   0.333333333333333   0.094117647058824   0.356862745098039
   0.337254901960784   0.094117647058824   0.356862745098039
   0.345098039215686   0.090196078431373   0.356862745098039
   0.349019607843137   0.094117647058824   0.360784313725490
   0.356862745098039   0.094117647058824   0.364705882352941
   0.360784313725490   0.094117647058824   0.364705882352941
   0.368627450980392   0.094117647058824   0.368627450980392
   0.372549019607843   0.094117647058824   0.368627450980392
   0.384313725490196   0.098039215686275   0.372549019607843
   0.384313725490196   0.098039215686275   0.372549019607843
   0.392156862745098   0.098039215686275   0.376470588235294
   0.403921568627451   0.101960784313725   0.380392156862745
   0.407843137254902   0.098039215686275   0.380392156862745
   0.415686274509804   0.101960784313725   0.380392156862745
   0.423529411764706   0.105882352941176   0.384313725490196
   0.435294117647059   0.109803921568627   0.388235294117647
   0.439215686274510   0.105882352941176   0.388235294117647
   0.447058823529412   0.109803921568627   0.392156862745098
   0.454901960784314   0.109803921568627   0.392156862745098
   0.466666666666667   0.117647058823529   0.396078431372549
   0.474509803921569   0.117647058823529   0.400000000000000
   0.478431372549020   0.117647058823529   0.400000000000000
   0.490196078431373   0.117647058823529   0.400000000000000
   0.501960784313725   0.125490196078431   0.403921568627451
   0.509803921568627   0.125490196078431   0.407843137254902
   0.521568627450980   0.129411764705882   0.407843137254902
   0.529411764705882   0.129411764705882   0.407843137254902
   0.541176470588235   0.137254901960784   0.415686274509804
   0.549019607843137   0.137254901960784   0.411764705882353
   0.552941176470588   0.141176470588235   0.411764705882353
   0.560784313725490   0.141176470588235   0.411764705882353
   0.564705882352941   0.141176470588235   0.407843137254902
   0.568627450980392   0.145098039215686   0.407843137254902
   0.576470588235294   0.145098039215686   0.403921568627451
   0.580392156862745   0.149019607843137   0.403921568627451
   0.584313725490196   0.149019607843137   0.403921568627451
   0.592156862745098   0.152941176470588   0.400000000000000
   0.596078431372549   0.152941176470588   0.400000000000000
   0.600000000000000   0.156862745098039   0.396078431372549
   0.607843137254902   0.156862745098039   0.396078431372549
   0.611764705882353   0.160784313725490   0.392156862745098
   0.615686274509804   0.160784313725490   0.392156862745098
   0.623529411764706   0.164705882352941   0.388235294117647
   0.627450980392157   0.164705882352941   0.384313725490196
   0.631372549019608   0.168627450980392   0.384313725490196
   0.639215686274510   0.172549019607843   0.380392156862745
   0.643137254901961   0.172549019607843   0.380392156862745
   0.647058823529412   0.176470588235294   0.376470588235294
   0.654901960784314   0.176470588235294   0.376470588235294
   0.658823529411765   0.180392156862745   0.372549019607843
   0.662745098039216   0.180392156862745   0.368627450980392
   0.670588235294118   0.184313725490196   0.368627450980392
   0.674509803921569   0.184313725490196   0.364705882352941
   0.678431372549020   0.188235294117647   0.360784313725490
   0.686274509803922   0.192156862745098   0.360784313725490
   0.690196078431373   0.192156862745098   0.356862745098039
   0.694117647058824   0.196078431372549   0.352941176470588
   0.701960784313725   0.196078431372549   0.352941176470588
   0.705882352941177   0.200000000000000   0.349019607843137
   0.709803921568627   0.203921568627451   0.345098039215686
   0.713725490196078   0.203921568627451   0.341176470588235
   0.721568627450980   0.207843137254902   0.341176470588235
   0.725490196078431   0.211764705882353   0.337254901960784
   0.729411764705882   0.211764705882353   0.333333333333333
   0.733333333333333   0.215686274509804   0.329411764705882
   0.741176470588235   0.219607843137255   0.325490196078431
   0.745098039215686   0.219607843137255   0.325490196078431
   0.749019607843137   0.223529411764706   0.321568627450980
   0.752941176470588   0.227450980392157   0.317647058823529
   0.760784313725490   0.231372549019608   0.313725490196078
   0.764705882352941   0.231372549019608   0.309803921568627
   0.768627450980392   0.235294117647059   0.309803921568627
   0.772549019607843   0.239215686274510   0.305882352941176
   0.776470588235294   0.243137254901961   0.301960784313725
   0.784313725490196   0.243137254901961   0.298039215686275
   0.788235294117647   0.247058823529412   0.294117647058824
   0.792156862745098   0.250980392156863   0.290196078431373
   0.796078431372549   0.254901960784314   0.286274509803922
   0.800000000000000   0.258823529411765   0.282352941176471
   0.803921568627451   0.262745098039216   0.278431372549020
   0.811764705882353   0.266666666666667   0.278431372549020
   0.815686274509804   0.266666666666667   0.274509803921569
   0.819607843137255   0.270588235294118   0.270588235294118
   0.823529411764706   0.274509803921569   0.266666666666667
   0.827450980392157   0.278431372549020   0.262745098039216
   0.831372549019608   0.282352941176471   0.258823529411765
   0.835294117647059   0.286274509803922   0.254901960784314
   0.839215686274510   0.290196078431373   0.250980392156863
   0.843137254901961   0.294117647058824   0.247058823529412
   0.847058823529412   0.298039215686275   0.243137254901961
   0.850980392156863   0.301960784313725   0.239215686274510
   0.854901960784314   0.305882352941176   0.235294117647059
   0.858823529411765   0.309803921568627   0.231372549019608
   0.862745098039216   0.313725490196078   0.227450980392157
   0.866666666666667   0.317647058823529   0.223529411764706
   0.870588235294118   0.325490196078431   0.219607843137255
   0.874509803921569   0.329411764705882   0.215686274509804
   0.878431372549020   0.333333333333333   0.211764705882353
   0.882352941176471   0.337254901960784   0.207843137254902
   0.886274509803922   0.341176470588235   0.203921568627451
   0.890196078431372   0.345098039215686   0.200000000000000
   0.894117647058824   0.349019607843137   0.196078431372549
   0.894117647058824   0.356862745098039   0.192156862745098
   0.898039215686275   0.360784313725490   0.188235294117647
   0.901960784313726   0.364705882352941   0.184313725490196
   0.905882352941176   0.368627450980392   0.180392156862745
   0.909803921568627   0.372549019607843   0.176470588235294
   0.909803921568627   0.380392156862745   0.172549019607843
   0.913725490196078   0.384313725490196   0.168627450980392
   0.917647058823529   0.388235294117647   0.164705882352941
   0.921568627450980   0.396078431372549   0.160784313725490
   0.921568627450980   0.400000000000000   0.156862745098039
   0.925490196078431   0.403921568627451   0.149019607843137
   0.929411764705882   0.407843137254902   0.145098039215686
   0.933333333333333   0.415686274509804   0.141176470588235
   0.933333333333333   0.419607843137255   0.137254901960784
   0.937254901960784   0.427450980392157   0.133333333333333
   0.937254901960784   0.431372549019608   0.129411764705882
   0.941176470588235   0.435294117647059   0.125490196078431
   0.945098039215686   0.443137254901961   0.121568627450980
   0.945098039215686   0.447058823529412   0.117647058823529
   0.949019607843137   0.450980392156863   0.113725490196078
   0.949019607843137   0.458823529411765   0.109803921568627
   0.952941176470588   0.462745098039216   0.101960784313725
   0.952941176470588   0.470588235294118   0.098039215686275
   0.956862745098039   0.474509803921569   0.094117647058824
   0.956862745098039   0.482352941176471   0.090196078431373
   0.960784313725490   0.486274509803922   0.086274509803922
   0.960784313725490   0.490196078431373   0.082352941176471
   0.964705882352941   0.498039215686275   0.078431372549020
   0.964705882352941   0.501960784313725   0.074509803921569
   0.968627450980392   0.509803921568627   0.066666666666667
   0.968627450980392   0.513725490196078   0.062745098039216
   0.968627450980392   0.521568627450980   0.058823529411765
   0.972549019607843   0.525490196078431   0.054901960784314
   0.972549019607843   0.533333333333333   0.050980392156863
   0.976470588235294   0.537254901960784   0.047058823529412
   0.976470588235294   0.545098039215686   0.043137254901961
   0.976470588235294   0.552941176470588   0.039215686274510
   0.976470588235294   0.556862745098039   0.035294117647059
   0.980392156862745   0.564705882352941   0.031372549019608
   0.980392156862745   0.568627450980392   0.031372549019608
   0.980392156862745   0.576470588235294   0.027450980392157
   0.980392156862745   0.580392156862745   0.027450980392157
   0.984313725490196   0.588235294117647   0.023529411764706
   0.984313725490196   0.592156862745098   0.023529411764706
   0.984313725490196   0.600000000000000   0.023529411764706
   0.984313725490196   0.607843137254902   0.023529411764706
   0.984313725490196   0.611764705882353   0.023529411764706
   0.984313725490196   0.619607843137255   0.027450980392157
   0.988235294117647   0.623529411764706   0.027450980392157
   0.988235294117647   0.631372549019608   0.031372549019608
   0.988235294117647   0.639215686274510   0.035294117647059
   0.988235294117647   0.643137254901961   0.039215686274510
   0.988235294117647   0.650980392156863   0.043137254901961
   0.988235294117647   0.654901960784314   0.047058823529412
   0.988235294117647   0.662745098039216   0.054901960784314
   0.988235294117647   0.670588235294118   0.058823529411765
   0.988235294117647   0.674509803921569   0.066666666666667
   0.988235294117647   0.682352941176471   0.070588235294118
   0.988235294117647   0.690196078431373   0.078431372549020
   0.988235294117647   0.694117647058824   0.086274509803922
   0.988235294117647   0.701960784313725   0.090196078431373
   0.988235294117647   0.709803921568627   0.098039215686275
   0.984313725490196   0.713725490196078   0.105882352941176
   0.984313725490196   0.721568627450980   0.113725490196078
   0.984313725490196   0.729411764705882   0.121568627450980
   0.984313725490196   0.733333333333333   0.129411764705882
   0.984313725490196   0.741176470588235   0.133333333333333
   0.984313725490196   0.749019607843137   0.141176470588235
   0.980392156862745   0.752941176470588   0.149019607843137
   0.980392156862745   0.760784313725490   0.156862745098039
   0.980392156862745   0.768627450980392   0.168627450980392
   0.980392156862745   0.772549019607843   0.176470588235294
   0.976470588235294   0.780392156862745   0.184313725490196
   0.976470588235294   0.788235294117647   0.192156862745098
   0.976470588235294   0.792156862745098   0.200000000000000
   0.972549019607843   0.800000000000000   0.207843137254902
   0.972549019607843   0.807843137254902   0.219607843137255
   0.972549019607843   0.811764705882353   0.227450980392157
   0.968627450980392   0.819607843137255   0.235294117647059
   0.968627450980392   0.827450980392157   0.247058823529412
   0.968627450980392   0.831372549019608   0.254901960784314
   0.964705882352941   0.839215686274510   0.266666666666667
   0.964705882352941   0.847058823529412   0.278431372549020
   0.960784313725490   0.850980392156863   0.286274509803922
   0.960784313725490   0.858823529411765   0.298039215686275
   0.960784313725490   0.866666666666667   0.309803921568627
   0.956862745098039   0.870588235294118   0.321568627450980
   0.956862745098039   0.878431372549020   0.329411764705882
   0.952941176470588   0.886274509803922   0.341176470588235
   0.952941176470588   0.890196078431372   0.352941176470588
   0.952941176470588   0.898039215686275   0.368627450980392
   0.949019607843137   0.901960784313726   0.380392156862745
   0.949019607843137   0.909803921568627   0.392156862745098
   0.949019607843137   0.913725490196078   0.403921568627451
   0.949019607843137   0.921568627450980   0.419607843137255
   0.945098039215686   0.925490196078431   0.431372549019608
   0.945098039215686   0.933333333333333   0.447058823529412
   0.945098039215686   0.937254901960784   0.458823529411765
   0.945098039215686   0.941176470588235   0.474509803921569
   0.949019607843137   0.949019607843137   0.486274509803922
   0.949019607843137   0.952941176470588   0.501960784313725
   0.949019607843137   0.956862745098039   0.513725490196078
   0.952941176470588   0.960784313725490   0.529411764705882
   0.956862745098039   0.968627450980392   0.545098039215686
   0.956862745098039   0.972549019607843   0.556862745098039
   0.960784313725490   0.976470588235294   0.568627450980392
   0.964705882352941   0.980392156862745   0.584313725490196
   0.968627450980392   0.984313725490196   0.596078431372549
   0.972549019607843   0.988235294117647   0.607843137254902
   0.976470588235294   0.992156862745098   0.619607843137255
   0.984313725490196   0.996078431372549   0.631372549019608
   0.988235294117647   1.000000000000000   0.643137254901961
    ];

P = size(values,1);
map = interp1(1:size(values,1), values, linspace(1,P,n), 'linear');
