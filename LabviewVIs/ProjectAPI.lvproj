<?xml version='1.0' encoding='UTF-8'?>
<Project Type="Project" LVVersion="19008000">
	<Item Name="My Computer" Type="My Computer">
		<Property Name="NI.SortType" Type="Int">3</Property>
		<Property Name="server.app.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.control.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.tcp.enabled" Type="Bool">false</Property>
		<Property Name="server.tcp.port" Type="Int">0</Property>
		<Property Name="server.tcp.serviceName" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.tcp.serviceName.default" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.vi.callsEnabled" Type="Bool">true</Property>
		<Property Name="server.vi.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="specify.custom.address" Type="Bool">false</Property>
		<Item Name="controls" Type="Folder">
			<Item Name="ActiveTracksFunctionalGlobalOperations.ctl" Type="VI" URL="../ActiveTracksFunctionalGlobalOperations.ctl"/>
			<Item Name="CalibrationStepsEnum.ctl" Type="VI" URL="../CalibrationStepsEnum.ctl"/>
			<Item Name="CircleDetectionEnum.ctl" Type="VI" URL="../CircleDetectionEnum.ctl"/>
			<Item Name="CompletedTracksFunctionalGlobalOperations.ctl" Type="VI" URL="../CompletedTracksFunctionalGlobalOperations.ctl"/>
			<Item Name="FauxWorm.ctl" Type="VI" URL="../FauxWorm.ctl"/>
			<Item Name="LinearTransformationDescriptionCluster.ctl" Type="VI" URL="../LinearTransformationDescriptionCluster.ctl"/>
			<Item Name="Parameters.ctl" Type="VI" URL="../Parameters.ctl"/>
			<Item Name="StimulusSnippetCluster.ctl" Type="VI" URL="../StimulusSnippetCluster.ctl"/>
			<Item Name="TimestampsFunctionalGlobalOperations.ctl" Type="VI" URL="../TimestampsFunctionalGlobalOperations.ctl"/>
			<Item Name="TrackTerminationReasonEnum.ctl" Type="VI" URL="../TrackTerminationReasonEnum.ctl"/>
			<Item Name="WormTrackCluster.ctl" Type="VI" URL="../WormTrackCluster.ctl"/>
			<Item Name="WormTrackSavingCluster.ctl" Type="VI" URL="../WormTrackSavingCluster.ctl"/>
			<Item Name="WormTrackSnippetCluster.ctl" Type="VI" URL="../WormTrackSnippetCluster.ctl"/>
		</Item>
		<Item Name="globals" Type="Folder">
			<Item Name="active_tracks_functional_global.vi" Type="VI" URL="../active_tracks_functional_global.vi"/>
			<Item Name="completed_tracks_functional_gloabl.vi" Type="VI" URL="../completed_tracks_functional_gloabl.vi"/>
			<Item Name="DebugGlobal.vi" Type="VI" URL="../DebugGlobal.vi"/>
			<Item Name="ParametersGlobal.vi" Type="VI" URL="../ParametersGlobal.vi"/>
			<Item Name="timestamps_functional_global.vi" Type="VI" URL="../timestamps_functional_global.vi"/>
		</Item>
		<Item Name="HDF5" Type="Folder">
			<Item Name="matlab_create.vi" Type="VI" URL="../matlab_create.vi"/>
			<Item Name="matlab_locale.vi" Type="VI" URL="../matlab_locale.vi"/>
			<Item Name="matlab_write.vi" Type="VI" URL="../matlab_write.vi"/>
			<Item Name="save_tracks_to_mat.vi" Type="VI" URL="../save_tracks_to_mat.vi"/>
		</Item>
		<Item Name="math_helper_functions" Type="Folder">
			<Item Name="1D_index_to_pixel_XY.vi" Type="VI" URL="../1D_index_to_pixel_XY.vi"/>
			<Item Name="absolute_orientation.vi" Type="VI" URL="../absolute_orientation.vi"/>
			<Item Name="affine_transform_point.vi" Type="VI" URL="../affine_transform_point.vi"/>
			<Item Name="array_mean_double.vi" Type="VI" URL="../array_mean_double.vi"/>
			<Item Name="array_median_double.vi" Type="VI" URL="../array_median_double.vi"/>
			<Item Name="calibration_image_to_conversion_ratio.vi" Type="VI" URL="../calibration_image_to_conversion_ratio.vi"/>
			<Item Name="CameraXY2projectorXY.vi" Type="VI" URL="../CameraXY2projectorXY.vi"/>
			<Item Name="center2rectangle.vi" Type="VI" URL="../center2rectangle.vi"/>
			<Item Name="circlefilter.vi" Type="VI" URL="../circlefilter.vi"/>
			<Item Name="convert_negative_blue_to_units_of_red.vi" Type="VI" URL="../convert_negative_blue_to_units_of_red.vi"/>
			<Item Name="enumerate_structuring_elements.vi" Type="VI" URL="../enumerate_structuring_elements.vi"/>
			<Item Name="euclidean_distance_double.vi" Type="VI" URL="../euclidean_distance_double.vi"/>
			<Item Name="euclidean_distance_single.vi" Type="VI" URL="../euclidean_distance_single.vi"/>
			<Item Name="FindClosestPoint.vi" Type="VI" URL="../FindClosestPoint.vi"/>
			<Item Name="FindClosestPointArray.vi" Type="VI" URL="../FindClosestPointArray.vi"/>
			<Item Name="flip_2D_single_array_td.vi" Type="VI" URL="../flip_2D_single_array_td.vi"/>
			<Item Name="interp_1d_nearest.vi" Type="VI" URL="../interp_1d_nearest.vi"/>
			<Item Name="linspace_double.vi" Type="VI" URL="../linspace_double.vi"/>
			<Item Name="mean_point_distance.vi" Type="VI" URL="../mean_point_distance.vi"/>
			<Item Name="normalize_vector_array_single.vi" Type="VI" URL="../normalize_vector_array_single.vi"/>
			<Item Name="parfor_progress.vi" Type="VI" URL="../parfor_progress.vi"/>
			<Item Name="pixel_XY_to_1D_index.vi" Type="VI" URL="../pixel_XY_to_1D_index.vi"/>
			<Item Name="points_to_semicolon_delimited_string.vi" Type="VI" URL="../points_to_semicolon_delimited_string.vi"/>
			<Item Name="ProjectorCameraXY2LinearTransformation.vi" Type="VI" URL="../ProjectorCameraXY2LinearTransformation.vi"/>
			<Item Name="projectorXY2CameraXY.vi" Type="VI" URL="../projectorXY2CameraXY.vi"/>
			<Item Name="overlay_crossmarks.vi" Type="VI" URL="../overlay_crossmarks.vi"/>
			<Item Name="radiucenter2boundingbox.vi" Type="VI" URL="../radiucenter2boundingbox.vi"/>
			<Item Name="radiucenter2topleftwidthheight.vi" Type="VI" URL="../radiucenter2topleftwidthheight.vi"/>
			<Item Name="radiucenterinbound.vi" Type="VI" URL="../radiucenterinbound.vi"/>
			<Item Name="random_array_element_single.vi" Type="VI" URL="../random_array_element_single.vi"/>
			<Item Name="ramp_by_delta_double.vi" Type="VI" URL="../ramp_by_delta_double.vi"/>
			<Item Name="reflect_and_rotate_point.vi" Type="VI" URL="../reflect_and_rotate_point.vi"/>
			<Item Name="rotate_and_reflect_point.vi" Type="VI" URL="../rotate_and_reflect_point.vi"/>
			<Item Name="rotate_point_by_radians_double.vi" Type="VI" URL="../rotate_point_by_radians_double.vi"/>
			<Item Name="rotate_point_by_radians_single.vi" Type="VI" URL="../rotate_point_by_radians_single.vi"/>
			<Item Name="sample_from_gaussian.vi" Type="VI" URL="../sample_from_gaussian.vi"/>
			<Item Name="semicolon_delimited_string_to_points.vi" Type="VI" URL="../semicolon_delimited_string_to_points.vi"/>
			<Item Name="sort_points_iterative_closest_point.vi" Type="VI" URL="../sort_points_iterative_closest_point.vi"/>
			<Item Name="SquareSpiral.vi" Type="VI" URL="../SquareSpiral.vi"/>
			<Item Name="stimulus_intensity_to_red_blue_color.vi" Type="VI" URL="../stimulus_intensity_to_red_blue_color.vi"/>
		</Item>
		<Item Name="centerlines" Type="Folder">
			<Item Name="iterative_segments_to_centerline.vi" Type="VI" URL="../iterative_segments_to_centerline.vi"/>
			<Item Name="iterative_segments_to_centerline_oneiteration.vi" Type="VI" URL="../iterative_segments_to_centerline_oneiteration.vi"/>
			<Item Name="iterative_thinned_image_to_segment.vi" Type="VI" URL="../iterative_thinned_image_to_segment.vi"/>
			<Item Name="resample_centerline.vi" Type="VI" URL="../resample_centerline.vi"/>
			<Item Name="select_largest_area_blob.vi" Type="VI" URL="../select_largest_area_blob.vi"/>
			<Item Name="thinned_image_to_centerline.vi" Type="VI" URL="../thinned_image_to_centerline.vi"/>
			<Item Name="thinned_image_to_centerline_debug.vi" Type="VI" URL="../thinned_image_to_centerline_debug.vi"/>
		</Item>
		<Item Name="calibration" Type="Folder">
			<Item Name="calibration_main.vi" Type="VI" URL="../calibration_main.vi"/>
			<Item Name="calculate_deviation_from_spatial_standard.vi" Type="VI" URL="../calculate_deviation_from_spatial_standard.vi"/>
			<Item Name="calculate_frames_per_second_count_dropped_frames.vi" Type="VI" URL="../calculate_frames_per_second_count_dropped_frames.vi"/>
			<Item Name="check_drift_during_experiment.vi" Type="VI" URL="../check_drift_during_experiment.vi"/>
			<Item Name="debug_temporal_calibration.vi" Type="VI" URL="../debug_temporal_calibration.vi"/>
			<Item Name="draw_faux_worm.vi" Type="VI" URL="../draw_faux_worm.vi"/>
			<Item Name="draw_temporal_calibration.vi" Type="VI" URL="../draw_temporal_calibration.vi"/>
			<Item Name="extract_calibration_standard_point_locations.vi" Type="VI" URL="../extract_calibration_standard_point_locations.vi"/>
			<Item Name="generate_spatial_calibration_pattern.vi" Type="VI" URL="../generate_spatial_calibration_pattern.vi"/>
			<Item Name="prompt_user_for_number_input.vi" Type="VI" URL="../prompt_user_for_number_input.vi"/>
			<Item Name="prompt_user_for_projector_settings.vi" Type="VI" URL="../prompt_user_for_projector_settings.vi"/>
			<Item Name="read_temporal_calibration.vi" Type="VI" URL="../read_temporal_calibration.vi"/>
			<Item Name="update_faux_worm_positions.vi" Type="VI" URL="../update_faux_worm_positions.vi"/>
			<Item Name="update_faux_worm_set.vi" Type="VI" URL="../update_faux_worm_set.vi"/>
			<Item Name="update_temporal_calibration_locations.vi" Type="VI" URL="../update_temporal_calibration_locations.vi"/>
			<Item Name="white_out_TC_zone.vi" Type="VI" URL="../white_out_TC_zone.vi"/>
			<Item Name="write_parameters_csv_from_calibration.vi" Type="VI" URL="../write_parameters_csv_from_calibration.vi"/>
		</Item>
		<Item Name="configuration" Type="Folder">
			<Item Name="display_projector_one_color.vi" Type="VI" URL="../display_projector_one_color.vi"/>
			<Item Name="get_camera_framerate.vi" Type="VI" URL="../get_camera_framerate.vi"/>
			<Item Name="get_git_hash.vi" Type="VI" URL="../get_git_hash.vi"/>
			<Item Name="initialize_projector.vi" Type="VI" URL="../initialize_projector.vi"/>
			<Item Name="load_parameters.vi" Type="VI" URL="../load_parameters.vi"/>
			<Item Name="load_parameters_csv.vi" Type="VI" URL="../load_parameters_csv.vi"/>
			<Item Name="projector_exposure_to_bit_depth.vi" Type="VI" URL="../projector_exposure_to_bit_depth.vi"/>
			<Item Name="setup_window_display.vi" Type="VI" URL="../setup_window_display.vi"/>
			<Item Name="snap_single_image.vi" Type="VI" URL="../snap_single_image.vi"/>
			<Item Name="test_connection.vi" Type="VI" URL="../test_connection.vi"/>
			<Item Name="turn_off_projector.vi" Type="VI" URL="../turn_off_projector.vi"/>
		</Item>
		<Item Name="post_processing" Type="Folder">
			<Item Name="encode_camera_frames.vi" Type="VI" URL="../encode_camera_frames.vi"/>
			<Item Name="get_all_experimental_folders_given_date.vi" Type="VI" URL="../get_all_experimental_folders_given_date.vi"/>
			<Item Name="projector_image_to_camera_image.vi" Type="VI" URL="../projector_image_to_camera_image.vi"/>
			<Item Name="projector_image_to_camera_image_subvi.vi" Type="VI" URL="../projector_image_to_camera_image_subvi.vi"/>
			<Item Name="post_processing.vi" Type="VI" URL="../post_processing.vi"/>
			<Item Name="transfer_and_submit_job.vi" Type="VI" URL="../transfer_and_submit_job.vi"/>
		</Item>
		<Item Name="tracking" Type="Folder">
			<Item Name="display_active_track_snippets.vi" Type="VI" URL="../display_active_track_snippets.vi"/>
			<Item Name="draw_stimulus_given_tracks.vi" Type="VI" URL="../draw_stimulus_given_tracks.vi"/>
			<Item Name="draw_tracking_image.vi" Type="VI" URL="../draw_tracking_image.vi"/>
			<Item Name="filter_active_tracks_for_mature_tracks.vi" Type="VI" URL="../filter_active_tracks_for_mature_tracks.vi"/>
			<Item Name="mature_track.vi" Type="VI" URL="../mature_track.vi"/>
			<Item Name="track2behavior.vi" Type="VI" URL="../track2behavior.vi"/>
			<Item Name="update_track_given_snippet.vi" Type="VI" URL="../update_track_given_snippet.vi"/>
			<Item Name="update_track_snippet.vi" Type="VI" URL="../update_track_snippet.vi"/>
		</Item>
		<Item Name="unused_helper_functions" Type="Folder">
			<Item Name="resave_camera_images_parfor.vi" Type="VI" URL="../resave_camera_images_parfor.vi"/>
			<Item Name="extract_calibration_standard_point_locations_debug.vi" Type="VI" URL="../extract_calibration_standard_point_locations_debug.vi"/>
		</Item>
		<Item Name="main_loop_functions" Type="Folder">
			<Item Name="main_image_aquisition.vi" Type="VI" URL="../main_image_aquisition.vi"/>
			<Item Name="main_save_camera_images.vi" Type="VI" URL="../main_save_camera_images.vi"/>
			<Item Name="main_median_filter_update.vi" Type="VI" URL="../main_median_filter_update.vi"/>
			<Item Name="main_draw_projector.vi" Type="VI" URL="../main_draw_projector.vi"/>
			<Item Name="main_save_projector_images.vi" Type="VI" URL="../main_save_projector_images.vi"/>
			<Item Name="main_track.vi" Type="VI" URL="../main_track.vi"/>
		</Item>
		<Item Name="HeadTailContinuousGWN.lvclass" Type="LVClass" URL="../HeadTailContinuousGWN.lvclass"/>
		<Item Name="main.vi" Type="VI" URL="../main.vi"/>
		<Item Name="StimulusObject.lvclass" Type="LVClass" URL="../StimulusObject.lvclass"/>
		<Item Name="parameters.csv" Type="Document" URL="../../parameters.csv"/>
		<Item Name="ConstantImageStimulus.lvclass" Type="LVClass" URL="../ConstantImageStimulusClass/ConstantImageStimulus.lvclass"/>
		<Item Name="FullWormRails.lvclass" Type="LVClass" URL="../FullWormRails.lvclass"/>
		<Item Name="FullWormRailsTriggeredBySmallEllipseRatio.lvclass" Type="LVClass" URL="../FullWormRailsTriggeredBySmallEllipseRatio.lvclass"/>
		<Item Name="RunConstantImageStimulus.vi" Type="VI" URL="../RunConstantImageStimulus.vi"/>
		<Item Name="RunFullWormRails.vi" Type="VI" URL="../RunFullWormRails.vi"/>
		<Item Name="RunHeadorTailContinuousGWN.vi" Type="VI" URL="../RunHeadorTailContinuousGWN.vi"/>
		<Item Name="RunHeadTailContinuousGWN.vi" Type="VI" URL="../RunHeadTailContinuousGWN.vi"/>
		<Item Name="RunRailsTriggeredByTurning.vi" Type="VI" URL="../RunRailsTriggeredByTurning.vi"/>
		<Item Name="RunHeadandTailRailswithDelays.vi" Type="VI" URL="../RunHeadandTailRailswithDelays.vi"/>
		<Item Name="RunHeadorTailRails.vi" Type="VI" URL="../RunHeadorTailRails.vi"/>
		<Item Name="RunJustTailRails.vi" Type="VI" URL="../RunJustTailRails.vi"/>
		<Item Name="HeadorTailRails.lvclass" Type="LVClass" URL="../HeadorTailRails.lvclass"/>
		<Item Name="JustTailRails.lvclass" Type="LVClass" URL="../JustTailRails.lvclass"/>
		<Item Name="HeadorTailContinuousGWN.lvclass" Type="LVClass" URL="../HeadorTailContinuousGWN.lvclass"/>
		<Item Name="HeadandTailRailswithDelay.lvclass" Type="LVClass" URL="../HeadandTailRailswithDelay.lvclass"/>
		<Item Name="FullWormGWN.lvclass" Type="LVClass" URL="../FullWormGWN.lvclass"/>
		<Item Name="RunFullWormContinuousGWN.vi" Type="VI" URL="../RunFullWormContinuousGWN.vi"/>
		<Item Name="Dependencies" Type="Dependencies">
			<Item Name="user.lib" Type="Folder">
				<Item Name="Array of VData to VArray__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Array of VData to VArray__ogtk.vi"/>
				<Item Name="Build Error Cluster__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/error/error.llb/Build Error Cluster__ogtk.vi"/>
				<Item Name="Cluster to Array of VData__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Cluster to Array of VData__ogtk.vi"/>
				<Item Name="Get Array Element TDEnum__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get Array Element TDEnum__ogtk.vi"/>
				<Item Name="Get Data Name from TD__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get Data Name from TD__ogtk.vi"/>
				<Item Name="Get Data Name__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get Data Name__ogtk.vi"/>
				<Item Name="Get Header from TD__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get Header from TD__ogtk.vi"/>
				<Item Name="Get Last PString__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get Last PString__ogtk.vi"/>
				<Item Name="Get PString__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get PString__ogtk.vi"/>
				<Item Name="Get TDEnum from Data__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get TDEnum from Data__ogtk.vi"/>
				<Item Name="Get Variant Attributes__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Get Variant Attributes__ogtk.vi"/>
				<Item Name="Parse String with TDs__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Parse String with TDs__ogtk.vi"/>
				<Item Name="Set Data Name__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Set Data Name__ogtk.vi"/>
				<Item Name="Split Cluster TD__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Split Cluster TD__ogtk.vi"/>
				<Item Name="Type Descriptor Enumeration__ogtk.ctl" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Type Descriptor Enumeration__ogtk.ctl"/>
				<Item Name="Type Descriptor Header__ogtk.ctl" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Type Descriptor Header__ogtk.ctl"/>
				<Item Name="Type Descriptor__ogtk.ctl" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Type Descriptor__ogtk.ctl"/>
				<Item Name="Variant to Header Info__ogtk.vi" Type="VI" URL="/&lt;userlib&gt;/_OpenG.lib/lvdata/lvdata.llb/Variant to Header Info__ogtk.vi"/>
			</Item>
			<Item Name="vi.lib" Type="Folder">
				<Item Name="Application Directory.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Application Directory.vi"/>
				<Item Name="BuildHelpPath.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/BuildHelpPath.vi"/>
				<Item Name="Calibration Reference Points.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Calibration.llb/Calibration Reference Points.ctl"/>
				<Item Name="Check if File or Folder Exists.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/libraryn.llb/Check if File or Folder Exists.vi"/>
				<Item Name="Check Special Tags.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Check Special Tags.vi"/>
				<Item Name="Clear Errors.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Clear Errors.vi"/>
				<Item Name="Close File+.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Close File+.vi"/>
				<Item Name="Color (U64)" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/Color (U64)"/>
				<Item Name="Color to RGB.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/colorconv.llb/Color to RGB.vi"/>
				<Item Name="compatReadText.vi" Type="VI" URL="/&lt;vilib&gt;/_oldvers/_oldvers.llb/compatReadText.vi"/>
				<Item Name="Convert property node font to graphics font.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Convert property node font to graphics font.vi"/>
				<Item Name="Correction Learn Setup.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Calibration.llb/Correction Learn Setup.ctl"/>
				<Item Name="Details Display Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Details Display Dialog.vi"/>
				<Item Name="DialogType.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/DialogType.ctl"/>
				<Item Name="DialogTypeEnum.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/DialogTypeEnum.ctl"/>
				<Item Name="Distortion Model Types.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Calibration.llb/Distortion Model Types.ctl"/>
				<Item Name="Distortion Model.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Calibration.llb/Distortion Model.ctl"/>
				<Item Name="Error Cluster From Error Code.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Error Cluster From Error Code.vi"/>
				<Item Name="Error Code Database.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Error Code Database.vi"/>
				<Item Name="Error Statistics.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Calibration.llb/Error Statistics.ctl"/>
				<Item Name="ErrWarn.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/ErrWarn.ctl"/>
				<Item Name="eventvkey.ctl" Type="VI" URL="/&lt;vilib&gt;/event_ctls.llb/eventvkey.ctl"/>
				<Item Name="Find First Error.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Find First Error.vi"/>
				<Item Name="Find Tag.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Find Tag.vi"/>
				<Item Name="Format Message String.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Format Message String.vi"/>
				<Item Name="General Error Handler Core CORE.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/General Error Handler Core CORE.vi"/>
				<Item Name="General Error Handler.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/General Error Handler.vi"/>
				<Item Name="Generate Temporary File Path.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/libraryn.llb/Generate Temporary File Path.vi"/>
				<Item Name="Get String Text Bounds.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Get String Text Bounds.vi"/>
				<Item Name="Get Text Rect.vi" Type="VI" URL="/&lt;vilib&gt;/picture/picture.llb/Get Text Rect.vi"/>
				<Item Name="GetHelpDir.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/GetHelpDir.vi"/>
				<Item Name="GetRTHostConnectedProp.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/GetRTHostConnectedProp.vi"/>
				<Item Name="H5D.lvlib" Type="Library" URL="/&lt;vilib&gt;/addons/h5labview2/dataset/H5D.lvlib"/>
				<Item Name="H5Equery.vi" Type="VI" URL="/&lt;vilib&gt;/addons/h5labview2/base/H5Equery.vi"/>
				<Item Name="H5F.lvlib" Type="Library" URL="/&lt;vilib&gt;/addons/h5labview2/file/H5F.lvlib"/>
				<Item Name="H5G.lvlib" Type="Library" URL="/&lt;vilib&gt;/addons/h5labview2/group/H5G.lvlib"/>
				<Item Name="H5Lexists.vi" Type="VI" URL="/&lt;vilib&gt;/addons/h5labview2/base/H5Lexists.vi"/>
				<Item Name="H5P.lvlib" Type="Library" URL="/&lt;vilib&gt;/addons/h5labview2/props/H5P.lvlib"/>
				<Item Name="hid_t.ctl" Type="VI" URL="/&lt;vilib&gt;/addons/h5labview2/base/hid_t.ctl"/>
				<Item Name="Image Type" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/Image Type"/>
				<Item Name="Image Unit" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/Image Unit"/>
				<Item Name="IMAQ Copy" Type="VI" URL="/&lt;vilib&gt;/vision/Management.llb/IMAQ Copy"/>
				<Item Name="IMAQ Create" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ Create"/>
				<Item Name="IMAQ Dispose" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ Dispose"/>
				<Item Name="IMAQ GetImageInfo" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ GetImageInfo"/>
				<Item Name="IMAQ GetImageSize" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ GetImageSize"/>
				<Item Name="IMAQ Image.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/IMAQ Image.ctl"/>
				<Item Name="IMAQ ImageToArray" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ ImageToArray"/>
				<Item Name="IMAQ Overlay Multiple Lines 2" Type="VI" URL="/&lt;vilib&gt;/vision/Overlay.llb/IMAQ Overlay Multiple Lines 2"/>
				<Item Name="IMAQ Overlay Oval" Type="VI" URL="/&lt;vilib&gt;/vision/Overlay.llb/IMAQ Overlay Oval"/>
				<Item Name="IMAQ Overlay Rectangle" Type="VI" URL="/&lt;vilib&gt;/vision/Overlay.llb/IMAQ Overlay Rectangle"/>
				<Item Name="IMAQ Overlay Text" Type="VI" URL="/&lt;vilib&gt;/vision/Overlay.llb/IMAQ Overlay Text"/>
				<Item Name="IMAQ Read Image And Vision Info 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files1.llb/IMAQ Read Image And Vision Info 2"/>
				<Item Name="IMAQ ReadFile 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ ReadFile 2"/>
				<Item Name="IMAQ SetImageSize" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ SetImageSize"/>
				<Item Name="IMAQ WindClose" Type="VI" URL="/&lt;vilib&gt;/vision/Display.llb/IMAQ WindClose"/>
				<Item Name="IMAQ WindDraw" Type="VI" URL="/&lt;vilib&gt;/vision/Display.llb/IMAQ WindDraw"/>
				<Item Name="IMAQ WindMove" Type="VI" URL="/&lt;vilib&gt;/vision/Display.llb/IMAQ WindMove"/>
				<Item Name="IMAQ WindSetup" Type="VI" URL="/&lt;vilib&gt;/vision/Display.llb/IMAQ WindSetup"/>
				<Item Name="IMAQ WindSize" Type="VI" URL="/&lt;vilib&gt;/vision/Display.llb/IMAQ WindSize"/>
				<Item Name="IMAQ WindZoom 2" Type="VI" URL="/&lt;vilib&gt;/vision/Display.llb/IMAQ WindZoom 2"/>
				<Item Name="IMAQ Write BMP File 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ Write BMP File 2"/>
				<Item Name="IMAQ Write File 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ Write File 2"/>
				<Item Name="IMAQ Write Image And Vision Info File 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ Write Image And Vision Info File 2"/>
				<Item Name="IMAQ Write JPEG File 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ Write JPEG File 2"/>
				<Item Name="IMAQ Write JPEG2000 File 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ Write JPEG2000 File 2"/>
				<Item Name="IMAQ Write PNG File 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ Write PNG File 2"/>
				<Item Name="IMAQ Write TIFF File 2" Type="VI" URL="/&lt;vilib&gt;/vision/Files.llb/IMAQ Write TIFF File 2"/>
				<Item Name="IMAQdx.ctl" Type="VI" URL="/&lt;vilib&gt;/userDefined/High Color/IMAQdx.ctl"/>
				<Item Name="lib_path.vi" Type="VI" URL="/&lt;vilib&gt;/addons/h5labview2/base/lib_path.vi"/>
				<Item Name="List Directory and LLBs.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/libraryn.llb/List Directory and LLBs.vi"/>
				<Item Name="Longest Line Length in Pixels.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Longest Line Length in Pixels.vi"/>
				<Item Name="LVBoundsTypeDef.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/miscctls.llb/LVBoundsTypeDef.ctl"/>
				<Item Name="LVRectTypeDef.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/miscctls.llb/LVRectTypeDef.ctl"/>
				<Item Name="NI_AALBase.lvlib" Type="Library" URL="/&lt;vilib&gt;/Analysis/NI_AALBase.lvlib"/>
				<Item Name="NI_FileType.lvlib" Type="Library" URL="/&lt;vilib&gt;/Utility/lvfile.llb/NI_FileType.lvlib"/>
				<Item Name="NI_Gmath.lvlib" Type="Library" URL="/&lt;vilib&gt;/gmath/NI_Gmath.lvlib"/>
				<Item Name="NI_PackedLibraryUtility.lvlib" Type="Library" URL="/&lt;vilib&gt;/Utility/LVLibp/NI_PackedLibraryUtility.lvlib"/>
				<Item Name="NI_PtbyPt.lvlib" Type="Library" URL="/&lt;vilib&gt;/ptbypt/NI_PtbyPt.lvlib"/>
				<Item Name="NI_Vision_Acquisition_Software.lvlib" Type="Library" URL="/&lt;vilib&gt;/vision/driver/NI_Vision_Acquisition_Software.lvlib"/>
				<Item Name="NI_Vision_Development_Module.lvlib" Type="Library" URL="/&lt;vilib&gt;/vision/NI_Vision_Development_Module.lvlib"/>
				<Item Name="Not Found Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Not Found Dialog.vi"/>
				<Item Name="Open File+.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Open File+.vi"/>
				<Item Name="Particle Parameters" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/Particle Parameters"/>
				<Item Name="Random Number (Range) DBL.vi" Type="VI" URL="/&lt;vilib&gt;/numeric/Random Number (Range) DBL.vi"/>
				<Item Name="Random Number (Range) I64.vi" Type="VI" URL="/&lt;vilib&gt;/numeric/Random Number (Range) I64.vi"/>
				<Item Name="Random Number (Range) U64.vi" Type="VI" URL="/&lt;vilib&gt;/numeric/Random Number (Range) U64.vi"/>
				<Item Name="Random Number (Range).vi" Type="VI" URL="/&lt;vilib&gt;/numeric/Random Number (Range).vi"/>
				<Item Name="Read Delimited Spreadsheet (DBL).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read Delimited Spreadsheet (DBL).vi"/>
				<Item Name="Read Delimited Spreadsheet (I64).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read Delimited Spreadsheet (I64).vi"/>
				<Item Name="Read Delimited Spreadsheet (string).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read Delimited Spreadsheet (string).vi"/>
				<Item Name="Read Delimited Spreadsheet.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read Delimited Spreadsheet.vi"/>
				<Item Name="Read File+ (string).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read File+ (string).vi"/>
				<Item Name="Read From Spreadsheet File (DBL).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read From Spreadsheet File (DBL).vi"/>
				<Item Name="Read From Spreadsheet File (I64).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read From Spreadsheet File (I64).vi"/>
				<Item Name="Read From Spreadsheet File (string).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read From Spreadsheet File (string).vi"/>
				<Item Name="Read From Spreadsheet File.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read From Spreadsheet File.vi"/>
				<Item Name="Read Lines From File (with error IO).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read Lines From File (with error IO).vi"/>
				<Item Name="Read Lines From File.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Read Lines From File.vi"/>
				<Item Name="Recursive File List.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/libraryn.llb/Recursive File List.vi"/>
				<Item Name="rel_path.vi" Type="VI" URL="/&lt;vilib&gt;/addons/h5labview2/base/rel_path.vi"/>
				<Item Name="Remove Duplicates From 1D Array.vim" Type="VI" URL="/&lt;vilib&gt;/Array/Remove Duplicates From 1D Array.vim"/>
				<Item Name="ROI Descriptor" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/ROI Descriptor"/>
				<Item Name="Search and Replace Pattern.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Search and Replace Pattern.vi"/>
				<Item Name="Set Bold Text.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Set Bold Text.vi"/>
				<Item Name="Set String Value.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Set String Value.vi"/>
				<Item Name="Simple Error Handler.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Simple Error Handler.vi"/>
				<Item Name="Simple Grid Descriptor" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/Simple Grid Descriptor"/>
				<Item Name="Space Constant.vi" Type="VI" URL="/&lt;vilib&gt;/dlg_ctls.llb/Space Constant.vi"/>
				<Item Name="sub_Random U32.vi" Type="VI" URL="/&lt;vilib&gt;/numeric/sub_Random U32.vi"/>
				<Item Name="System Exec.vi" Type="VI" URL="/&lt;vilib&gt;/Platform/system.llb/System Exec.vi"/>
				<Item Name="TagReturnType.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/TagReturnType.ctl"/>
				<Item Name="Three Button Dialog CORE.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Three Button Dialog CORE.vi"/>
				<Item Name="Three Button Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Three Button Dialog.vi"/>
				<Item Name="Trim Whitespace.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/Trim Whitespace.vi"/>
				<Item Name="whitespace.ctl" Type="VI" URL="/&lt;vilib&gt;/Utility/error.llb/whitespace.ctl"/>
				<Item Name="Write Delimited Spreadsheet (DBL).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Write Delimited Spreadsheet (DBL).vi"/>
				<Item Name="Write Delimited Spreadsheet (I64).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Write Delimited Spreadsheet (I64).vi"/>
				<Item Name="Write Delimited Spreadsheet (string).vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Write Delimited Spreadsheet (string).vi"/>
				<Item Name="Write Delimited Spreadsheet.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Write Delimited Spreadsheet.vi"/>
				<Item Name="Write Spreadsheet String.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/file.llb/Write Spreadsheet String.vi"/>
				<Item Name="Path To Command Line String.vi" Type="VI" URL="/&lt;vilib&gt;/AdvancedString/Path To Command Line String.vi"/>
				<Item Name="New Zip File.vi" Type="VI" URL="/&lt;vilib&gt;/zip/New Zip File.vi"/>
				<Item Name="Relative Path To Platform Independent String.vi" Type="VI" URL="/&lt;vilib&gt;/AdvancedString/Relative Path To Platform Independent String.vi"/>
				<Item Name="Add File to Zip.vi" Type="VI" URL="/&lt;vilib&gt;/zip/Add File to Zip.vi"/>
				<Item Name="Compare Two Paths.vi" Type="VI" URL="/&lt;vilib&gt;/Utility/libraryn.llb/Compare Two Paths.vi"/>
				<Item Name="Close Zip File.vi" Type="VI" URL="/&lt;vilib&gt;/zip/Close Zip File.vi"/>
				<Item Name="PathToUNIXPathString.vi" Type="VI" URL="/&lt;vilib&gt;/Platform/CFURL.llb/PathToUNIXPathString.vi"/>
				<Item Name="IMAQ ColorImageToArray" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ ColorImageToArray"/>
				<Item Name="Delimited String to 1D String Array.vi" Type="VI" URL="/&lt;vilib&gt;/AdvancedString/Delimited String to 1D String Array.vi"/>
				<Item Name="ex_CorrectErrorChain.vi" Type="VI" URL="/&lt;vilib&gt;/express/express shared/ex_CorrectErrorChain.vi"/>
				<Item Name="subDisplayMessage.vi" Type="VI" URL="/&lt;vilib&gt;/express/express output/DisplayMessageBlock.llb/subDisplayMessage.vi"/>
				<Item Name="IMAQ ArrayToColorImage" Type="VI" URL="/&lt;vilib&gt;/vision/Basics.llb/IMAQ ArrayToColorImage"/>
				<Item Name="subFile Dialog.vi" Type="VI" URL="/&lt;vilib&gt;/express/express input/FileDialogBlock.llb/subFile Dialog.vi"/>
				<Item Name="1D String Array to Delimited String.vi" Type="VI" URL="/&lt;vilib&gt;/AdvancedString/1D String Array to Delimited String.vi"/>
				<Item Name="IMAQ Rounding Mode.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/IMAQ Rounding Mode.ctl"/>
				<Item Name="Vision Info Type2.ctl" Type="VI" URL="/&lt;vilib&gt;/vision/Image Controls.llb/Vision Info Type2.ctl"/>
				<Item Name="IMAQ Overlay Line" Type="VI" URL="/&lt;vilib&gt;/vision/Overlay.llb/IMAQ Overlay Line"/>
				<Item Name="IMAQ Copy Vision Info" Type="VI" URL="/&lt;vilib&gt;/vision/Management.llb/IMAQ Copy Vision Info"/>
			</Item>
			<Item Name="niimaqdx.dll" Type="Document" URL="niimaqdx.dll">
				<Property Name="NI.PreserveRelativePath" Type="Bool">true</Property>
			</Item>
			<Item Name="nivision.dll" Type="Document" URL="nivision.dll">
				<Property Name="NI.PreserveRelativePath" Type="Bool">true</Property>
			</Item>
			<Item Name="nivissvc.dll" Type="Document" URL="nivissvc.dll">
				<Property Name="NI.PreserveRelativePath" Type="Bool">true</Property>
			</Item>
			<Item Name="random_array_element_interger.vi" Type="VI" URL="../random_array_element_interger.vi"/>
		</Item>
		<Item Name="Build Specifications" Type="Build">
			<Item Name="projector_image_to_camera_image_worker" Type="EXE">
				<Property Name="App_copyErrors" Type="Bool">true</Property>
				<Property Name="App_INI_aliasGUID" Type="Str">{DD838B18-FA1B-4735-8D66-96BF93402484}</Property>
				<Property Name="App_INI_GUID" Type="Str">{B4362BF2-90A6-485D-8A96-4A25305EE6B3}</Property>
				<Property Name="App_serverConfig.httpPort" Type="Int">8002</Property>
				<Property Name="Bld_autoIncrement" Type="Bool">true</Property>
				<Property Name="Bld_buildCacheID" Type="Str">{518825EB-502C-45A7-BEEB-11E79D98763D}</Property>
				<Property Name="Bld_buildSpecName" Type="Str">projector_image_to_camera_image_worker</Property>
				<Property Name="Bld_excludeInlineSubVIs" Type="Bool">true</Property>
				<Property Name="Bld_excludeLibraryItems" Type="Bool">true</Property>
				<Property Name="Bld_excludePolymorphicVIs" Type="Bool">true</Property>
				<Property Name="Bld_localDestDir" Type="Path">../builds/NI_AB_PROJECTNAME/projector_image_to_camera_image_worker</Property>
				<Property Name="Bld_localDestDirType" Type="Str">relativeToCommon</Property>
				<Property Name="Bld_modifyLibraryFile" Type="Bool">true</Property>
				<Property Name="Bld_previewCacheID" Type="Str">{944EF78F-B678-4EDC-80BE-B8D0BC755772}</Property>
				<Property Name="Bld_version.build" Type="Int">2</Property>
				<Property Name="Bld_version.major" Type="Int">1</Property>
				<Property Name="Destination[0].destName" Type="Str">projector_image_to_camera_image_worker.exe</Property>
				<Property Name="Destination[0].path" Type="Path">../builds/NI_AB_PROJECTNAME/projector_image_to_camera_image_worker/projector_image_to_camera_image_worker.exe</Property>
				<Property Name="Destination[0].preserveHierarchy" Type="Bool">true</Property>
				<Property Name="Destination[0].type" Type="Str">App</Property>
				<Property Name="Destination[1].destName" Type="Str">Support Directory</Property>
				<Property Name="Destination[1].path" Type="Path">../builds/NI_AB_PROJECTNAME/projector_image_to_camera_image_worker/data</Property>
				<Property Name="DestinationCount" Type="Int">2</Property>
				<Property Name="Exe_cmdLineArgs" Type="Bool">true</Property>
				<Property Name="Source[0].itemID" Type="Str">{5AF22950-F1EB-4A9C-8C7B-BA12267398B9}</Property>
				<Property Name="Source[0].type" Type="Str">Container</Property>
				<Property Name="Source[1].destinationIndex" Type="Int">0</Property>
				<Property Name="Source[1].itemID" Type="Ref"></Property>
				<Property Name="Source[1].sourceInclusion" Type="Str">TopLevel</Property>
				<Property Name="Source[1].type" Type="Str">VI</Property>
				<Property Name="SourceCount" Type="Int">2</Property>
				<Property Name="TgtF_companyName" Type="Str">Princeton University</Property>
				<Property Name="TgtF_fileDescription" Type="Str">projector_image_to_camera_image_worker</Property>
				<Property Name="TgtF_internalName" Type="Str">projector_image_to_camera_image_worker</Property>
				<Property Name="TgtF_legalCopyright" Type="Str">Copyright © 2020 Princeton University</Property>
				<Property Name="TgtF_productName" Type="Str">projector_image_to_camera_image_worker</Property>
				<Property Name="TgtF_targetfileGUID" Type="Str">{B60C3873-C414-4AE0-AA2F-4E8ABD16EC5B}</Property>
				<Property Name="TgtF_targetfileName" Type="Str">projector_image_to_camera_image_worker.exe</Property>
				<Property Name="TgtF_versionIndependent" Type="Bool">true</Property>
			</Item>
		</Item>
	</Item>
</Project>
