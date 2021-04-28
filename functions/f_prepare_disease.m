% NOT IMPLEMENTED YET !!!!!

function [POP] = prepare_disease(POP, CORAL, META)

%     % MODIFY DISEASE PREVALENCE IN THE PARAMETER FILE
%     % numcorals1=find(pop(:,:,6:10)>0); %pulls out all corals
%     % numcorals2=length(numcorals1); % use this one if prevalence of all corals

%     numcorals8a=find(pop(:,:,8)>0); %pulls out all corals
%     numcorals8b=length(numcorals8a);
%     numcorals9a=find(pop(:,:,9)>0); %pulls out all corals
%     numcorals9b=length(numcorals9a);
%     numcorals10a=find(pop(:,:,10)>0); %pulls out all corals
%     numcorals10b=length(numcorals10a);
%     
%     %BBD % use code 11 to indicate was diseased but died
%     num_bbd_8_1 = find(pop(:,:,11)==1 |pop(:,:,11)==11);% first work out prevalence for each vul coral type
%     num_bbd_8_2 = length(num_bbd_8_1);
%     current_bbd_prev_8=num_bbd_8_2/numcorals8b;
%     num_bbd_9_1 = find(pop(:,:,12)==1 | pop(:,:,12)==11);
%     num_bbd_9_2 = length(num_bbd_9_1);
%     current_bbd_prev_9=num_bbd_9_2/numcorals9b;
%     num_bbd_10_1 = find(pop(:,:,13)==1 | pop(:,:,13)==11);
%     num_bbd_10_2 = length(num_bbd_10_1);
%     current_bbd_prev_10=num_bbd_10_2/numcorals10b;
%     current_bbd_prev=[current_bbd_prev_8,current_bbd_prev_9,current_bbd_prev_10];
%     revised_bbd = zeros(1,3); % new bbd probabilities
%     for i=1:3
%         if BBD_PREV(1,i)>current_bbd_prev(1,i) % if inadequate disease prevalence
%             revised_bbd(1,i)=BBD_PREV(1,i)-current_bbd_prev(1,i);
%         end %therefore only raise disease if too low, if too high leave it alone (natural and would occur if non bbd corals died)
%     end
% 
%     %YBD
%     num_ybd_8_1 = find(pop(:,:,11)==2| pop(:,:,11)==12);% first work out prevalence for each vul coral type
%     num_ybd_8_2 = length(num_ybd_8_1);
%     current_ybd_prev_8=num_ybd_8_2/numcorals8b;
%     num_ybd_9_1 = find(pop(:,:,12)==2| pop(:,:,12)==12);
%     num_ybd_9_2 = length(num_ybd_9_1);
%     current_ybd_prev_9=num_ybd_9_2/numcorals9b;
%     num_ybd_10_1 = find(pop(:,:,13)==2 | pop(:,:,13)==12);
%     num_ybd_10_2 = length(num_ybd_10_1);
%     current_ybd_prev_10=num_ybd_10_2/numcorals10b;
%     current_ybd_prev=[current_ybd_prev_8,current_ybd_prev_9,current_ybd_prev_10];
%     revised_ybd = zeros(1,3); % new ybd probabilities
%     for i=1:3
%         if YBD_PREV(1,i)>current_ybd_prev(1,i) % if inadequate disease prevalence
%             revised_ybd(1,i)=YBD_PREV(1,i)-current_ybd_prev(1,i);
%         end %therefore only raise disease if too low, if too high leave it alone (natural and would occur if non ybd corals died)
%     end
% 
%     %WPII
%     num_wpii_8_1 = find(pop(:,:,11)==3 | pop(:,:,11)==13);% first work out prevalence for each vul coral type
%     num_wpii_8_2 = length(num_wpii_8_1);
%     current_wpii_prev_8=num_wpii_8_2/numcorals8b;
%     num_wpii_9_1 = find(pop(:,:,12)==3 | pop(:,:,12)==13);
%     num_wpii_9_2 = length(num_wpii_9_1);
%     current_wpii_prev_9=num_wpii_9_2/numcorals9b;
%     num_wpii_10_1 = find(pop(:,:,13)==3 | pop(:,:,13)==13);
%     num_wpii_10_2 = length(num_wpii_10_1);
%     current_wpii_prev_10=num_wpii_10_2/numcorals10b;
%     current_wpii_prev=[current_wpii_prev_8,current_wpii_prev_9,current_wpii_prev_10];
%     revised_wpii = zeros(1,3); % new wpii probabilities
%     for i=1:3
%         if WPII_PREV(1,i)>current_wpii_prev(1,i) % if inadequate disease prevalence
%             revised_wpii(1,i)=WPII_PREV(1,i)-current_wpii_prev(1,i);
%         end %therefore only raise disease if too low, if too high leave it alone (natural and would occur if non bbd corals died)
%     end
% 
%     %WBD
%     num_wbd_8_1 = find(pop(:,:,11)==4 |pop(:,:,11)==14);% first work out prevalence for each vul coral type
%     num_wbd_8_2 = length(num_wbd_8_1);
%     current_wbd_prev_8=num_wbd_8_2/numcorals8b;
%     num_wbd_9_1 = find(pop(:,:,12)==4 | pop(:,:,12)==14);
%     num_wbd_9_2 = length(num_wbd_9_1);
%     current_wbd_prev_9=num_wbd_9_2/numcorals9b;
%     num_wbd_10_1 = find(pop(:,:,13)==4 | pop(:,:,13)==14);
%     num_wbd_10_2 = length(num_wbd_10_1);
%     current_wbd_prev_10=num_wbd_10_2/numcorals10b;
%     current_wbd_prev=[current_wbd_prev_8,current_wbd_prev_9,current_wbd_prev_10];
%     revised_wbd = zeros(1,3); % new wbd probabilities
%     for i=1:3
%         if WBD_PREV(1,i)>current_wbd_prev(1,i) % if inadequate disease prevalence
%             revised_wbd(1,i)=WBD_PREV(1,i)-current_wbd_prev(1,i);
%         end %therefore only raise disease if too low, if too high leave it alone (natural and would occur if non wbd corals died)
%     end
%     dis=[current_bbd_prev(:,:);current_ybd_prev(:,:);current_wpii_prev(:,:);...
%         current_wbd_prev(:,:)]; % monitor disease prevalence
% 
%     disease=[revised_bbd(:,:);revised_ybd(:,:);revised_wpii(:,:);...
%         revised_wbd(:,:)]; % send down to process cell
% end