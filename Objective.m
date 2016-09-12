function f=Objective(x)
f1=-0.013409341616*x(0+1)-0.013424916342*x(1+1)-0.012214099922*x(2+1)-0.007888735922*x(3+1)-0.013393766889*x(4+1)-0.014651307489*x(5+1)-0.015328647501*x(6+1)-0.011452107843*x(7+1)-0.012120651563*x(8+1)-0.016457899224*x(9+1)-0.022997074752*x(10+1)-0.022591047949*x(11+1)-0.006510288575*x(12+1)-0.016062563092*x(13+1)-0.037610704332*x(14+1)-0.055915009202*x(15+1)+0.008652348929*x(16+1)-0.003671360236*x(17+1)-0.055468170282*x(18+1)-0.163458284529*x(19+1)+0.036781699714*x(20+1)+0.064268046383*x(21+1)+0.143618696789*x(22+1)+0.261305185445*x(23+1)+0.036861127298*x(24+1)+0.064357904856*x(25+1)+0.143729297770*x(26+1)+0.261415129097*x(27+1)+0.008720868956*x(28+1)-0.003600870084*x(29+1)-0.055398136365*x(30+1)-0.163405888501*x(31+1)-0.006454646231*x(32+1)-0.016009014904*x(33+1)-0.037564055843*x(34+1)-0.055885402732*x(35+1)-0.012075792745*x(36+1)-0.016416487457*x(37+1)-0.022963669371*x(38+1)-0.022571666583*x(39+1)-0.013356244547*x(40+1)-0.014617472808*x(41+1)-0.015302467600*x(42+1)-0.011437594231*x(43+1)-0.013375468087*x(44+1)-0.013394691628*x(45+1)-0.012191133989*x(46+1)-0.007876242740*x(47+1)+0.000000082964;
% f2=-0.013409341616*x(96+1)-0.013424916342*x(97+1)-0.012214099922*x(98+1)-0.007888735922*x(99+1)-0.013393766889*x(100+1)-0.014651307489*x(101+1)-0.015328647501*x(102+1)-0.011452107843*x(103+1)-0.012120651563*x(104+1)-0.016457899224*x(105+1)-0.022997074752*x(106+1)-0.022591047949*x(107+1)-0.006510288575*x(108+1)-0.016062563092*x(109+1)-0.037610704332*x(110+1)-0.055915009202*x(111+1)+0.008652348929*x(112+1)-0.003671360236*x(113+1)-0.055468170282*x(114+1)-0.163458284529*x(115+1)+0.036781699714*x(116+1)+0.064268046383*x(117+1)+0.143618696789*x(118+1)+0.261305185445*x(119+1)+0.036861127298*x(120+1)+0.064357904856*x(121+1)+0.143729297770*x(122+1)+0.261415129097*x(123+1)+0.008720868956*x(124+1)-0.003600870084*x(125+1)-0.055398136365*x(126+1)-0.163405888501*x(127+1)-0.006454646231*x(128+1)-0.016009014904*x(129+1)-0.037564055843*x(130+1)-0.055885402732*x(131+1)-0.012075792745*x(132+1)-0.016416487457*x(133+1)-0.022963669371*x(134+1)-0.022571666583*x(135+1)-0.013356244547*x(136+1)-0.014617472808*x(137+1)-0.015302467600*x(138+1)-0.011437594231*x(139+1)-0.013375468087*x(140+1)-0.013394691628*x(141+1)-0.012191133989*x(142+1)-0.007876242740*x(143+1)-0.011387616314;
% f3=-0.013409341616*x(192+1)-0.013424916342*x(193+1)-0.012214099922*x(194+1)-0.007888735922*x(195+1)-0.013393766889*x(196+1)-0.014651307489*x(197+1)-0.015328647501*x(198+1)-0.011452107843*x(199+1)-0.012120651563*x(200+1)-0.016457899224*x(201+1)-0.022997074752*x(202+1)-0.022591047949*x(203+1)-0.006510288575*x(204+1)-0.016062563092*x(205+1)-0.037610704332*x(206+1)-0.055915009202*x(207+1)+0.008652348929*x(208+1)-0.003671360236*x(209+1)-0.055468170282*x(210+1)-0.163458284529*x(211+1)+0.036781699714*x(212+1)+0.064268046383*x(213+1)+0.143618696789*x(214+1)+0.261305185445*x(215+1)+0.036861127298*x(216+1)+0.064357904856*x(217+1)+0.143729297770*x(218+1)+0.261415129097*x(219+1)+0.008720868956*x(220+1)-0.003600870084*x(221+1)-0.055398136365*x(222+1)-0.163405888501*x(223+1)-0.006454646231*x(224+1)-0.016009014904*x(225+1)-0.037564055843*x(226+1)-0.055885402732*x(227+1)-0.012075792745*x(228+1)-0.016416487457*x(229+1)-0.022963669371*x(230+1)-0.022571666583*x(231+1)-0.013356244547*x(232+1)-0.014617472808*x(233+1)-0.015302467600*x(234+1)-0.011437594231*x(235+1)-0.013375468087*x(236+1)-0.013394691628*x(237+1)-0.012191133989*x(238+1)-0.007876242740*x(239+1)-0.008460848900;
f4=-0.013409341616*x(24+1)-0.013424916342*x(25+1)-0.012214099922*x(26+1)-0.007888735922*x(27+1)-0.013393766889*x(28+1)-0.014651307489*x(29+1)-0.015328647501*x(30+1)-0.011452107843*x(31+1)-0.012120651563*x(32+1)-0.016457899224*x(33+1)-0.022997074752*x(34+1)-0.022591047949*x(35+1)-0.006510288575*x(36+1)-0.016062563092*x(37+1)-0.037610704332*x(38+1)-0.055915009202*x(39+1)+0.008652348929*x(40+1)-0.003671360236*x(41+1)-0.055468170282*x(42+1)-0.163458284529*x(43+1)+0.036781699714*x(44+1)+0.064268046383*x(45+1)+0.143618696789*x(46+1)+0.261305185445*x(47+1)+0.036861127298*x(48+1)+0.064357904856*x(49+1)+0.143729297770*x(50+1)+0.261415129097*x(51+1)+0.008720868956*x(52+1)-0.003600870084*x(53+1)-0.055398136365*x(54+1)-0.163405888501*x(55+1)-0.006454646231*x(56+1)-0.016009014904*x(57+1)-0.037564055843*x(58+1)-0.055885402732*x(59+1)-0.012075792745*x(60+1)-0.016416487457*x(61+1)-0.022963669371*x(62+1)-0.022571666583*x(63+1)-0.013356244547*x(64+1)-0.014617472808*x(65+1)-0.015302467600*x(66+1)-0.011437594231*x(67+1)-0.013375468087*x(68+1)-0.013394691628*x(69+1)-0.012191133989*x(70+1)-0.007876242740*x(71+1)-0.008460848900;
% f5=-0.013409341616*x(120+1)-0.013424916342*x(121+1)-0.012214099922*x(122+1)-0.007888735922*x(123+1)-0.013393766889*x(124+1)-0.014651307489*x(125+1)-0.015328647501*x(126+1)-0.011452107843*x(127+1)-0.012120651563*x(128+1)-0.016457899224*x(129+1)-0.022997074752*x(130+1)-0.022591047949*x(131+1)-0.006510288575*x(132+1)-0.016062563092*x(133+1)-0.037610704332*x(134+1)-0.055915009202*x(135+1)+0.008652348929*x(136+1)-0.003671360236*x(137+1)-0.055468170282*x(138+1)-0.163458284529*x(139+1)+0.036781699714*x(140+1)+0.064268046383*x(141+1)+0.143618696789*x(142+1)+0.261305185445*x(143+1)+0.036861127298*x(144+1)+0.064357904856*x(145+1)+0.143729297770*x(146+1)+0.261415129097*x(147+1)+0.008720868956*x(148+1)-0.003600870084*x(149+1)-0.055398136365*x(150+1)-0.163405888501*x(151+1)-0.006454646231*x(152+1)-0.016009014904*x(153+1)-0.037564055843*x(154+1)-0.055885402732*x(155+1)-0.012075792745*x(156+1)-0.016416487457*x(157+1)-0.022963669371*x(158+1)-0.022571666583*x(159+1)-0.013356244547*x(160+1)-0.014617472808*x(161+1)-0.015302467600*x(162+1)-0.011437594231*x(163+1)-0.013375468087*x(164+1)-0.013394691628*x(165+1)-0.012191133989*x(166+1)-0.007876242740*x(167+1)-0.011387616314;
% f6=-0.013409341616*x(216+1)-0.013424916342*x(217+1)-0.012214099922*x(218+1)-0.007888735922*x(219+1)-0.013393766889*x(220+1)-0.014651307489*x(221+1)-0.015328647501*x(222+1)-0.011452107843*x(223+1)-0.012120651563*x(224+1)-0.016457899224*x(225+1)-0.022997074752*x(226+1)-0.022591047949*x(227+1)-0.006510288575*x(228+1)-0.016062563092*x(229+1)-0.037610704332*x(230+1)-0.055915009202*x(231+1)+0.008652348929*x(232+1)-0.003671360236*x(233+1)-0.055468170282*x(234+1)-0.163458284529*x(235+1)+0.036781699714*x(236+1)+0.064268046383*x(237+1)+0.143618696789*x(238+1)+0.261305185445*x(239+1)+0.036861127298*x(240+1)+0.064357904856*x(241+1)+0.143729297770*x(242+1)+0.261415129097*x(243+1)+0.008720868956*x(244+1)-0.003600870084*x(245+1)-0.055398136365*x(246+1)-0.163405888501*x(247+1)-0.006454646231*x(248+1)-0.016009014904*x(249+1)-0.037564055843*x(250+1)-0.055885402732*x(251+1)-0.012075792745*x(252+1)-0.016416487457*x(253+1)-0.022963669371*x(254+1)-0.022571666583*x(255+1)-0.013356244547*x(256+1)-0.014617472808*x(257+1)-0.015302467600*x(258+1)-0.011437594231*x(259+1)-0.013375468087*x(260+1)-0.013394691628*x(261+1)-0.012191133989*x(262+1)-0.007876242740*x(263+1)-0.000000082964;
f7=-0.274163010139*x(48+1)-0.244630511357*x(49+1)-0.185879931058*x(50+1)-0.101116369547*x(51+1)-0.303695508920*x(52+1)-0.273848592875*x(53+1)-0.211892912269*x(54+1)-0.117469177585*x(55+1)-0.363074923747*x(56+1)-0.335175438953*x(57+1)-0.270373947558*x(58+1)-0.156867428522*x(59+1)-0.450353823368*x(60+1)-0.433404291634*x(61+1)-0.377560010487*x(62+1)-0.239626588946*x(63+1)-0.554582254723*x(64+1)-0.570527893726*x(65+1)-0.566835213812*x(66+1)-0.424078916773*x(67+1)-0.642865047075*x(68+1)-0.727289814735*x(69+1)-0.895174034262*x(70+1)-0.889853864335*x(71+1)+0.642865047056*x(72+1)+0.727289814713*x(73+1)+0.895174034238*x(74+1)+0.889853864326*x(75+1)+0.554582254707*x(76+1)+0.570527893710*x(77+1)+0.566835213798*x(78+1)+0.424078916766*x(79+1)+0.450353823356*x(80+1)+0.433404291622*x(81+1)+0.377560010478*x(82+1)+0.239626588941*x(83+1)+0.363074923737*x(84+1)+0.335175438945*x(85+1)+0.270373947551*x(86+1)+0.156867428518*x(87+1)+0.303695508912*x(88+1)+0.273848592868*x(89+1)+0.211892912263*x(90+1)+0.117469177582*x(91+1)+0.274163010131*x(92+1)+0.244630511351*x(93+1)+0.185879931053*x(94+1)+0.101116369545*x(95+1)-0.001342978395;
% f8=-0.274163010139*x(144+1)-0.244630511357*x(145+1)-0.185879931058*x(146+1)-0.101116369547*x(147+1)-0.303695508920*x(148+1)-0.273848592875*x(149+1)-0.211892912269*x(150+1)-0.117469177585*x(151+1)-0.363074923747*x(152+1)-0.335175438953*x(153+1)-0.270373947558*x(154+1)-0.156867428522*x(155+1)-0.450353823368*x(156+1)-0.433404291634*x(157+1)-0.377560010487*x(158+1)-0.239626588946*x(159+1)-0.554582254723*x(160+1)-0.570527893726*x(161+1)-0.566835213812*x(162+1)-0.424078916773*x(163+1)-0.642865047075*x(164+1)-0.727289814735*x(165+1)-0.895174034262*x(166+1)-0.889853864335*x(167+1)+0.642865047056*x(168+1)+0.727289814713*x(169+1)+0.895174034238*x(170+1)+0.889853864326*x(171+1)+0.554582254707*x(172+1)+0.570527893710*x(173+1)+0.566835213798*x(174+1)+0.424078916766*x(175+1)+0.450353823356*x(176+1)+0.433404291622*x(177+1)+0.377560010478*x(178+1)+0.239626588941*x(179+1)+0.363074923737*x(180+1)+0.335175438945*x(181+1)+0.270373947551*x(182+1)+0.156867428518*x(183+1)+0.303695508912*x(184+1)+0.273848592868*x(185+1)+0.211892912263*x(186+1)+0.117469177582*x(187+1)+0.274163010131*x(188+1)+0.244630511351*x(189+1)+0.185879931053*x(190+1)+0.101116369545*x(191+1)+0.000000003889;
% f9=-0.274163010139*x(240+1)-0.244630511357*x(241+1)-0.185879931058*x(242+1)-0.101116369547*x(243+1)-0.303695508920*x(244+1)-0.273848592875*x(245+1)-0.211892912269*x(246+1)-0.117469177585*x(247+1)-0.363074923747*x(248+1)-0.335175438953*x(249+1)-0.270373947558*x(250+1)-0.156867428522*x(251+1)-0.450353823368*x(252+1)-0.433404291634*x(253+1)-0.377560010487*x(254+1)-0.239626588946*x(255+1)-0.554582254723*x(256+1)-0.570527893726*x(257+1)-0.566835213812*x(258+1)-0.424078916773*x(259+1)-0.642865047075*x(260+1)-0.727289814735*x(261+1)-0.895174034262*x(262+1)-0.889853864335*x(263+1)+0.642865047056*x(264+1)+0.727289814713*x(265+1)+0.895174034238*x(266+1)+0.889853864326*x(267+1)+0.554582254707*x(268+1)+0.570527893710*x(269+1)+0.566835213798*x(270+1)+0.424078916766*x(271+1)+0.450353823356*x(272+1)+0.433404291622*x(273+1)+0.377560010478*x(274+1)+0.239626588941*x(275+1)+0.363074923737*x(276+1)+0.335175438945*x(277+1)+0.270373947551*x(278+1)+0.156867428518*x(279+1)+0.303695508912*x(280+1)+0.273848592868*x(281+1)+0.211892912263*x(282+1)+0.117469177582*x(283+1)+0.274163010131*x(284+1)+0.244630511351*x(285+1)+0.185879931053*x(286+1)+0.101116369545*x(287+1)-0.000000002765;
f10=-0.013409341616*x(48+1)-0.013424916342*x(49+1)-0.012214099922*x(50+1)-0.007888735922*x(51+1)-0.013393766889*x(52+1)-0.014651307489*x(53+1)-0.015328647501*x(54+1)-0.011452107843*x(55+1)-0.012120651563*x(56+1)-0.016457899224*x(57+1)-0.022997074752*x(58+1)-0.022591047949*x(59+1)-0.006510288575*x(60+1)-0.016062563092*x(61+1)-0.037610704332*x(62+1)-0.055915009202*x(63+1)+0.008652348929*x(64+1)-0.003671360236*x(65+1)-0.055468170282*x(66+1)-0.163458284529*x(67+1)+0.036781699714*x(68+1)+0.064268046383*x(69+1)+0.143618696789*x(70+1)+0.261305185445*x(71+1)+0.036861127298*x(72+1)+0.064357904856*x(73+1)+0.143729297770*x(74+1)+0.261415129097*x(75+1)+0.008720868956*x(76+1)-0.003600870084*x(77+1)-0.055398136365*x(78+1)-0.163405888501*x(79+1)-0.006454646231*x(80+1)-0.016009014904*x(81+1)-0.037564055843*x(82+1)-0.055885402732*x(83+1)-0.012075792745*x(84+1)-0.016416487457*x(85+1)-0.022963669371*x(86+1)-0.022571666583*x(87+1)-0.013356244547*x(88+1)-0.014617472808*x(89+1)-0.015302467600*x(90+1)-0.011437594231*x(91+1)-0.013375468087*x(92+1)-0.013394691628*x(93+1)-0.012191133989*x(94+1)-0.007876242740*x(95+1)-0.000000082964;
% f11=-0.013409341616*x(144+1)-0.013424916342*x(145+1)-0.012214099922*x(146+1)-0.007888735922*x(147+1)-0.013393766889*x(148+1)-0.014651307489*x(149+1)-0.015328647501*x(150+1)-0.011452107843*x(151+1)-0.012120651563*x(152+1)-0.016457899224*x(153+1)-0.022997074752*x(154+1)-0.022591047949*x(155+1)-0.006510288575*x(156+1)-0.016062563092*x(157+1)-0.037610704332*x(158+1)-0.055915009202*x(159+1)+0.008652348929*x(160+1)-0.003671360236*x(161+1)-0.055468170282*x(162+1)-0.163458284529*x(163+1)+0.036781699714*x(164+1)+0.064268046383*x(165+1)+0.143618696789*x(166+1)+0.261305185445*x(167+1)+0.036861127298*x(168+1)+0.064357904856*x(169+1)+0.143729297770*x(170+1)+0.261415129097*x(171+1)+0.008720868956*x(172+1)-0.003600870084*x(173+1)-0.055398136365*x(174+1)-0.163405888501*x(175+1)-0.006454646231*x(176+1)-0.016009014904*x(177+1)-0.037564055843*x(178+1)-0.055885402732*x(179+1)-0.012075792745*x(180+1)-0.016416487457*x(181+1)-0.022963669371*x(182+1)-0.022571666583*x(183+1)-0.013356244547*x(184+1)-0.014617472808*x(185+1)-0.015302467600*x(186+1)-0.011437594231*x(187+1)-0.013375468087*x(188+1)-0.013394691628*x(189+1)-0.012191133989*x(190+1)-0.007876242740*x(191+1)-0.011387616314;
% f12=-0.013409341616*x(240+1)-0.013424916342*x(241+1)-0.012214099922*x(242+1)-0.007888735922*x(243+1)-0.013393766889*x(244+1)-0.014651307489*x(245+1)-0.015328647501*x(246+1)-0.011452107843*x(247+1)-0.012120651563*x(248+1)-0.016457899224*x(249+1)-0.022997074752*x(250+1)-0.022591047949*x(251+1)-0.006510288575*x(252+1)-0.016062563092*x(253+1)-0.037610704332*x(254+1)-0.055915009202*x(255+1)+0.008652348929*x(256+1)-0.003671360236*x(257+1)-0.055468170282*x(258+1)-0.163458284529*x(259+1)+0.036781699714*x(260+1)+0.064268046383*x(261+1)+0.143618696789*x(262+1)+0.261305185445*x(263+1)+0.036861127298*x(264+1)+0.064357904856*x(265+1)+0.143729297770*x(266+1)+0.261415129097*x(267+1)+0.008720868956*x(268+1)-0.003600870084*x(269+1)-0.055398136365*x(270+1)-0.163405888501*x(271+1)-0.006454646231*x(272+1)-0.016009014904*x(273+1)-0.037564055843*x(274+1)-0.055885402732*x(275+1)-0.012075792745*x(276+1)-0.016416487457*x(277+1)-0.022963669371*x(278+1)-0.022571666583*x(279+1)-0.013356244547*x(280+1)-0.014617472808*x(281+1)-0.015302467600*x(282+1)-0.011437594231*x(283+1)-0.013375468087*x(284+1)-0.013394691628*x(285+1)-0.012191133989*x(286+1)-0.007876242740*x(287+1)+0.008460848900;
f=f1+f4+f7+f10;
end