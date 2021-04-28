% NOT IMPLEMENTED YET !!!!!

% Dans POP(x,y).disease(s,c) ajouter des 0, 1, 2, 3 ou 4 pour les record

function POPxy = f_coral_disease(POPxy, s, META, REEF, CORAL)
%             if (c==8 & pop(x,y,c)>=CORAL3_MIN_SIZE)|(c==9 & pop(x,y,c)>=CORAL4_MIN_SIZE)|(c==10 & pop(x,y,c)>=CORAL5_MIN_SIZE)
%             % only vulnerable coral types of appropriate size
%             dis=pop(x,y,c+3);%checks corresponding layer 8-10 where value 0 none, 1 bbd, 2 ybd, 3 wpii, 4 wbd
%             if dis>9 % if it still has a code of 'died last time' - i.e., 11,12,13,14
%                 dis=0;
%             end
%             if dis==0 % if no disease then try to give all four disease
%                 for d=1:4 % move down 4 x 3 disease prob matrix of disease type x c
%                     if rand(1)<disease(d,c-7) %give it disease
%                         pop(x,y,c+6)=START_LESION_SIZE; % start disease recording
%                         pop(x,y,c)=pop(x,y,c)-START_LESION_SIZE; % kill coral
%                         pop(x,y,1)=pop(x,y,1)+START_LESION_SIZE; % add to turf
%                         pop(x,y,c+3)=d; % record that colony diseased
%                         break; % don't try to give it other diseases
%                     end
%                 end
%             elseif dis==1 % if bbd% already has disease now need to grow lesion
%                     newlesion = round(pi*((sqrt(pop(x,y,c+6)./pi))+BBD_PROG)^2); % new lesion size
%                     if (pop(x,y,c)-newlesion)<=0 % killed coral completely
%                         pop(x,y,1)=pop(x,y,1)+pop(x,y,c); % convert to turf
%                         pop(x,y,c)=0; % kill coral
%                         pop(x,y,5)=pop(x,y,5)-1;% remove record of coral
%                         pop(x,y,c+3)=11;% indicate died
%                         pop(x,y,c+6)=0;
%                         pop(x,y,c+11)=0; % remove bleaching history
%                     else
%                         pop(x,y,c)=pop(x,y,c)-newlesion;% remove tissue
%                         pop(x,y,1)=pop(x,y,1)+newlesion;% add to turf
%                         pop(x,y,c+6)=newlesion; % record new state
%                     end
%              elseif dis==2 % YBD
%                         newlesion = round(pi*((sqrt(pop(x,y,c+6)./pi))+YBD_PROG)^2); % new lesion size
%                     if (pop(x,y,c)-newlesion)<=0 % killed coral completely
%                         pop(x,y,1)=pop(x,y,1)+pop(x,y,c); % convert to turf
%                         pop(x,y,c)=0; % kill coral
%                         pop(x,y,5)=pop(x,y,5)-1;% remove record of coral
%                         pop(x,y,c+11)=0; % remove bleaching history
%                         pop(x,y,c+3)=12;% remove disease recording
%                         pop(x,y,c+6)=0;
%                     else
%                         pop(x,y,c)=pop(x,y,c)-newlesion;% remove tissue
%                         pop(x,y,1)=pop(x,y,1)+newlesion;% add to turf
%                         pop(x,y,c+6)=newlesion; % record new state
%                     end
%               elseif dis==3 % WPII
%                        newlesion = round(pi*((sqrt(pop(x,y,c+6)./pi))+WPII_PROG)^2); % new lesion size
%                     if (pop(x,y,c)-newlesion)<=0 % killed coral completely
%                         pop(x,y,1)=pop(x,y,1)+pop(x,y,c); % convert to turf
%                         pop(x,y,c)=0; % kill coral
%                         pop(x,y,5)=pop(x,y,5)-1;% remove record of coral
%                         pop(x,y,c+11)=0; % remove bleaching history
%                         pop(x,y,c+3)=13;% remove disease recording
%                         pop(x,y,c+6)=0;
%                     else
%                         pop(x,y,c)=pop(x,y,c)-newlesion;% remove tissue
%                         pop(x,y,1)=pop(x,y,1)+newlesion;% add to turf
%                         pop(x,y,c+6)=newlesion; % record new state
%                     end
%               elseif dis==4
%                     newlesion = round(pi*((sqrt(pop(x,y,c+6)./pi))+WBD_PROG)^2); % new lesion size
%                     if (pop(x,y,c)-newlesion)<=0 % killed coral completely
%                         pop(x,y,1)=pop(x,y,1)+pop(x,y,c); % convert to turf
%                         pop(x,y,c)=0; % kill coral
%                         pop(x,y,5)=pop(x,y,5)-1;% remove record of coral
%                         pop(x,y,c+11)=0; % remove bleaching history
%                         pop(x,y,c+3)=14;% remove disease recording
%                         pop(x,y,c+6)=0;
%                     else
%                         pop(x,y,c)=pop(x,y,c)-newlesion;% remove tissue
%                         pop(x,y,1)=pop(x,y,1)+newlesion;% add to turf
%                         pop(x,y,c+6)=newlesion; % record new state
%                     end
%               %  end
%             end % end the elseif for dis==0
%             end % end the if vulnerable coral and size
%         end % end disease